/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "lss.hpp"

#define LSS_BLOCKSIZE 512
#define LSS_DGDG_TYPE 0
#define LSS_FVFV_TYPE 1
#define LSS_DGFV_TYPE 2

lss_t& lss_t::Setup(mesh_t& mesh, linAlg_t& linAlg,
                                 lssSettings_t& settings){

  lss_t* lss = new lss_t(mesh, linAlg, settings);

  lss->cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;
  lss->advection  = (settings.compareSetting("ADVECTION SOLVER", "TRUE")) ? 1:0;
  lss->redistance = (settings.compareSetting("REDISTANCE SOLVER", "TRUE")) ? 1:0;

  lss->subcellStabilization = (settings.compareSetting("STABILIZATION", "SUBCELL")) ? 1:0;

  //setup cubature
  if (lss->cubature) {
    mesh.CubatureSetup();
    mesh.CubatureNodes();
  }

  lss->NVfields = mesh.dim; 
  lss->Nfields  = (lss->redistance) ? 2:1; 
  dlong Nlocal  = mesh.Nelements*mesh.Np;
  dlong Nhalo   = mesh.totalHaloPairs*mesh.Np;

  //setup linear algebra module
  lss->linAlg.InitKernels({"innerProd","axpy","amx"}, mesh.comm);

  /*setup trace halo exchange */
  lss->traceHalo = mesh.HaloTraceSetup(lss->Nfields);

  // compute samples of q at interpolation nodes
  lss->q = (dfloat*) calloc((Nlocal+Nhalo)*lss->Nfields, sizeof(dfloat));
  lss->o_q = mesh.device.malloc((Nlocal+Nhalo)*lss->Nfields*sizeof(dfloat), lss->q); 


  // this holds the gradient of qM, qP
  lss->gradq = (dfloat*) calloc(Nlocal*lss->Nfields*mesh.dim, sizeof(dfloat));
  lss->o_gradq = mesh.device.malloc(Nlocal*lss->Nfields*mesh.dim*sizeof(dfloat),lss->gradq);
  
  // this saves the history for qP, and qM to second order recontruction
  if(settings.compareSetting("TIME RECONSTRUCTION", "ENO2")){
  lss->Nrecon = 4; // Second order reconstruction, 
  }else if(settings.compareSetting("TIME RECONSTRUCTION", "ENO3")){
  lss->Nrecon = 6; 
  }
  // AK: this could be there be three lets keep it simple for now
  lss->phiH = (dfloat *) calloc(Nlocal*lss->Nfields*lss->Nrecon,sizeof(dfloat));
  lss->o_phiH = mesh.device.malloc(Nlocal*lss->Nfields*lss->Nrecon*sizeof(dfloat),lss->phiH);

  // No that the actual solution is always a scalar!!!
  lss->phi  = (dfloat *) calloc(Nlocal, sizeof(dfloat)); 
  lss->o_phi  = mesh.device.malloc(Nlocal*sizeof(dfloat), lss->phi);
 

  lss->U   = (dfloat*) calloc((Nlocal+Nhalo)*lss->NVfields, sizeof(dfloat));
  lss->o_U = mesh.device.malloc((Nlocal+Nhalo)*lss->NVfields*sizeof(dfloat), lss->U);
  
  //storage for M*q during reporting
  lss->o_Mq = mesh.device.malloc((Nlocal+Nhalo)*lss->Nfields*sizeof(dfloat), lss->q);
  
  // OCCA build stuff
  occa::properties &kernelInfo = lss->props; 

  //add boundary data to kernel info
  string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;
  kernelInfo["defines/" "p_blockSize"]= (int)LSS_BLOCKSIZE; 

  kernelInfo["defines/" "p_Nfields"] = lss->Nfields;
  kernelInfo["defines/" "p_NVfields"]= lss->NVfields;
  kernelInfo["defines/" "p_Nrecon"]= lss->Nrecon;

  int maxNodes = mymax(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int NblockV = 1024/mesh.Np; // works for CUDA
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo["defines/" "p_NblockS"]= NblockS;
  
  if(lss->cubature){
    int cubNblockV = 512/mesh.cubNp; // works for CUDA
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubMaxNodes = mymax(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;

    int cubMaxNodes1 = mymax(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockS = 512/cubMaxNodes; // works for CUDA
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }

  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";


  // set kernel name suffix
  char *suffix;
  if(mesh.elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(mesh.elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(mesh.elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(mesh.elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];
  

  // kernels from lss advection volume file
  if(lss->cubature){
    sprintf(fileName, DLSS "/okl/lssAdvectionCubVolume%s.okl", suffix);
    sprintf(kernelName, "lssAdvectionCubatureVolume%s", suffix);
  }else{
    sprintf(fileName, DLSS "/okl/lssAdvectionVolume%s.okl", suffix);
    sprintf(kernelName, "lssAdvectionVolume%s", suffix);
  }

  lss->advectionVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);
  // kernels from surface file
  if(lss->cubature){
  sprintf(fileName, DLSS "/okl/lssAdvectionCubSurface%s.okl", suffix);
  sprintf(kernelName, "lssAdvectionCubatureSurface%s", suffix);
  }else{
  sprintf(fileName, DLSS "/okl/lssAdvectionSurface%s.okl", suffix);
  sprintf(kernelName, "lssAdvectionSurface%s", suffix);
  }

  lss->advectionSurfaceKernel = buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);

  

  if (mesh.dim==2) {
    sprintf(fileName, DLSS "/okl/lssInitialCondition2D.okl");
    if(lss->advection){
      sprintf(kernelName, "lssAdvectionInitialCondition2D");
    }
    else{
      sprintf(kernelName, "lssRedistanceInitialCondition2D");
    }
  } else {
    sprintf(fileName, DLSS "/okl/lssInitialCondition3D.okl");
    sprintf(kernelName, "lssInitialCondition3D");
  }

  lss->initialConditionKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);
  if(lss->redistance){
  sprintf(fileName, DLSS "/okl/lssInitialCondition2D.okl");
  sprintf(kernelName, "lssRedistanceSetAuxiliaryField2D");
  lss->setAuxiliaryFieldKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);
  }
  
  // This takes place of flow solver for simple problems
  if(lss->advection){
    if (mesh.dim==2) {
    sprintf(fileName, DLSS "/okl/lssSetFlowField2D.okl");
    sprintf(kernelName, "lssSetFlowField2D");
  } else {
    sprintf(fileName, DLSS "/okl/lssSetFlowField3D.okl");
    sprintf(kernelName, "lssSetFlowField3D");
  }
  lss->setFlowFieldKernel = buildKernel(mesh.device, fileName, kernelName,
                                                  kernelInfo, mesh.comm);
  }

   sprintf(fileName, DLSS "/okl/lssReconstructENO.okl");

   if(settings.compareSetting("TIME RECONSTRUCTION", "ENO2")){ 
    sprintf(kernelName, "lssReconstructENO2");
    lss->reconstructENOKernel = buildKernel(mesh.device, fileName, kernelName,
                                                kernelInfo, mesh.comm);
    sprintf(kernelName, "lssInitialHistoryENO2");
    lss->initialHistoryKernel = buildKernel(mesh.device, fileName, kernelName,
                                                kernelInfo, mesh.comm);
   }else{
    sprintf(kernelName, "lssReconstructENO3");
    lss->reconstructENOKernel = buildKernel(mesh.device, fileName, kernelName,
                                                kernelInfo, mesh.comm);
    sprintf(kernelName, "lssInitialHistoryENO3");
    lss->initialHistoryKernel = buildKernel(mesh.device, fileName, kernelName,
                                                kernelInfo, mesh.comm);
  }
// mass matrix operator
  // sprintf(fileName, LIBP_DIR "/core/okl/MassMatrixOperator%s.okl", suffix);
  sprintf(kernelName, "MassMatrixOperator%s", suffix);

  lss->MassMatrixKernel = buildKernel(mesh.device, fileName, kernelName,
                                            kernelInfo, mesh.comm);


   sprintf(fileName, DLSS "/okl/lssRegularizedSign2D.okl");
   sprintf(kernelName, "lssRegularizedSign2D");
   lss->regularizedSignKernel = buildKernel(mesh.device, fileName, kernelName,
                                                kernelInfo, mesh.comm);
   
  sprintf(fileName, DLSS "/okl/lssRedistanceVolume%s.okl", suffix);
  sprintf(kernelName, "lssRedistanceVolume%s", suffix);

  lss->redistanceVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);


  sprintf(fileName, DLSS "/okl/lssRedistanceSurface%s.okl", suffix);
  sprintf(kernelName, "lssRedistanceSurface%s", suffix);

  lss->redistanceSurfaceKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);


  sprintf(fileName, DLSS "/okl/lssGradient%s.okl", suffix);
  sprintf(kernelName, "lssGradientVolume%s", suffix);

  lss->gradientVolumeKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);

  sprintf(kernelName, "lssGradientSurface%s", suffix);
  lss->gradientSurfaceKernel =  buildKernel(mesh.device, fileName, kernelName,
                                         kernelInfo, mesh.comm);



// 
  if(lss->subcellStabilization){
    lss->SetupStabilizer();   

#if 1
    const dlong sNlocal = mesh.Nelements*lss->subcell->Np*lss->Nfields; 

  lss->svq   = (dfloat*) calloc(sNlocal, sizeof(dfloat));
  lss->o_svq = mesh.device.malloc(sNlocal*sizeof(dfloat), lss->svq); 

  lss->gsq = (dfloat*) calloc(Nlocal, sizeof(dfloat));
  lss->o_gsq = mesh.device.malloc(Nlocal*sizeof(dfloat), lss->gsq); 
#endif

  }


  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    lss->timeStepper = new TimeStepper::ab3(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *lss);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    if(lss->subcellStabilization)
    lss->timeStepper = new TimeStepper::lserk4_subcell(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, lss->Nfields, lss->subcell->Nsubcells, lss->Nfields, *lss);
    else
    lss->timeStepper = new TimeStepper::lserk4(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, lss->Nfields, *lss);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    lss->timeStepper = new TimeStepper::dopri5(mesh.Nelements, mesh.totalHaloPairs,
                                              mesh.Np, 1, *lss);
  }

  // // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat cfl = 0.25; // depends on the stability region size

  dfloat dt = cfl*hmin/((mesh.N+1.)*(mesh.N+1.));
  lss->timeStepper->SetTimeStep(dt);  
  lss->eps    = 1.f*hmin;  


  #if 1

   dfloat *degree = (dfloat*) calloc(Nlocal, sizeof(dfloat));

   for(int i=0; i<Nlocal; i++ ){
    degree[i] = 1.0; 
   }
   lss->o_invDegree = mesh.device.malloc(Nlocal*sizeof(dfloat),degree);

   mesh.ogs->GatherScatter(lss->o_invDegree, ogs_dfloat, ogs_add, ogs_sym); 

   lss->o_invDegree.copyTo(degree);

   for(int i=0; i<Nlocal; i++ ){
    degree[i] = 1./degree[i]; 
    // printf("%.4f\n", degree[i]);
   }

   lss->o_invDegree.copyFrom(degree);


  #endif

  return *lss;
}

lss_t::~lss_t() {
  advectionVolumeKernel.free();
  advectionSurfaceKernel.free();
  redistanceVolumeKernel.free();
  redistanceSurfaceKernel.free();

  MassMatrixKernel.free();
  initialConditionKernel.free();
  setFlowFieldKernel.free(); 
  regularizedSignKernel.free(); 

  redistanceVolumeKernel.free();
  redistanceSurfaceKernel.free();

  if(subcellStabilization){
    indicatorMDAKernel.free();
    findNeighKernel.free();
    projectKernel.free();
    reconstructKernel.free();
    partialRedistanceVolumeKernel.free();
    reconstructInternalFaceKernel.free();
    projectDGKernel.free();
    reconstructExternalFaceKernel.free();
    partialRedistanceSurfaceKernel.free();
    mixedRedistanceSurfaceKernel.free();
    subcellComputeKernel.free();
    if(subcell) delete subcell;
  }
  
  if (timeStepper) delete timeStepper;
  if (traceHalo) traceHalo->Free();
}


