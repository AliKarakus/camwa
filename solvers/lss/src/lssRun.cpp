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

void lss_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  if(advection){
  initialConditionKernel(mesh.Nelements,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);

  setFlowFieldKernel(mesh.Nelements,
                    startTime,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    o_U, 
                    o_q);
}else{

shiftIndex = 0; 
historyIndex = 0; 

rtime = (dfloat *)calloc(Nrecon, sizeof(dfloat)); 
o_rtime = device.malloc(Nrecon*sizeof(dfloat), rtime); 

rtime[shiftIndex] = startTime; 

initialConditionKernel(mesh.Nelements,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_phi);

setAuxiliaryFieldKernel(mesh.Nelements,
           startTime,
           mesh.o_x,
           mesh.o_y,
           mesh.o_z,
           o_phi,
           o_q);

 const dlong N = mesh.Nelements*mesh.Np*Nfields; 
 o_phiH.copyFrom(o_q, N*sizeof(dfloat),shiftIndex*N*sizeof(dfloat),0);
 shiftIndex = (shiftIndex+Nrecon-1)%Nrecon;


DetectTroubledCells(o_q, subcell->o_ElementList); 

// for(int i=0; i<Nrecon; i++){
// // Create History
// const dlong N = mesh.Nelements*mesh.Np*Nfields; 
// o_phiH.copyFrom(o_q, N*sizeof(dfloat), shiftIndex*N); 

// timeStepper->Run(o_q, i*idt, (i+1)*idt);

// }





 //  // Compute gradient for regularization
 //  traceHalo->ExchangeStart(o_q, 1, ogs_dfloat);

 //   gradientVolumeKernel(mesh.Nelements,
 //                       offset, 
 //                       mesh.o_vgeo,
 //                       mesh.o_DWmatrices, // o_Dmatrices
 //                       o_q,
 //                       o_gradq); 

 //  traceHalo->ExchangeFinish(o_q, 1, ogs_dfloat);

 //  gradientSurfaceKernel(mesh.Nelements,
 //                        offset, 
 //                        mesh.o_sgeo,
 //                        mesh.o_LIFTT,
 //                        mesh.o_vmapM,
 //                        mesh.o_vmapP,
 //                        mesh.o_EToB,
 //                        0,
 //                        mesh.o_x,
 //                        mesh.o_y,
 //                        mesh.o_z,
 //                        o_q,
 //                        o_gradq);


 //  regularizedSignKernel(mesh.Nelements, 
 //                       startTime,
 //                       eps, 
 //                       mesh.o_x,
 //                       mesh.o_y,
 //                       mesh.o_z,
 //                       o_q,
 //                       o_gradq, 
 //                       o_sgnq); 

 // if(subcellStabilization){


 //   const int all = 1; // all elements, do not use element list
 //    projectKernel(mesh.Nelements,
 //                  all,  
 //                  subcell->o_ElementList, 
 //                  subcell->o_PMT,
 //                  o_sgnq, 
 //                  o_ssgnq);


 //   DetectTroubledCells(o_q, subcell->o_ElementList); 

  
 //  }

}

 // dfloat outputInterval; 
 // settings.getSetting("OUTPUT INTERVAL", outputInterval);
// printf("%.4e\n",timeStepper->GetTimeStep() );
// finalTime = startTime + 12*timeStepper->GetTimeStep();
// printf("%.4e %.4e\n",finalTime, timeStepper->GetTimeStep() );



 // finalTime = startTime + 53*outputInterval;
timeStepper->Run(o_q, startTime, finalTime);

#if 0
shiftIndex = 0; 
historyIndex = 0; 

rtime = (dfloat *)calloc(Nrecon, sizeof(dfloat)); 
o_rtime = device.malloc(Nrecon*sizeof(dfloat), rtime); 

rtime[shiftIndex] = startTime; 

setAuxiliaryFieldKernel(mesh.Nelements,
           startTime,
           mesh.o_x,
           mesh.o_y,
           mesh.o_z,
           o_phi,
           o_q);

const dlong N = mesh.Nelements*mesh.Np*Nfields; 
 o_phiH.copyFrom(o_q, N*sizeof(dfloat),shiftIndex*N*sizeof(dfloat),0);
 shiftIndex = (shiftIndex+Nrecon-1)%Nrecon;



// finalTime = startTime + 53*outputInterval;
timeStepper->Run(o_q, startTime, finalTime);
#endif



}
