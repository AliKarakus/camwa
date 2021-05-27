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

#ifndef LSS_HPP
#define LSS_HPP 1

#include "core.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"
#include "subcell.hpp"

#define DLSS LIBP_DIR"/solvers/lss/"

class lssSettings_t: public settings_t {
public:
  lssSettings_t(MPI_Comm& _comm);
  void report();
  void parseFromFile(occaSettings_t& occaSettings,
                     meshSettings_t& meshSettings,
                     const string filename);
};



class lss_t: public solver_t {
public:
  int cubature; 
  int advection;
  int redistance;
  int subcellStabilization;
  int Nfields; 
  int NVfields; 
  int Nrecon; // required number of history point for reconstruction in time
  dlong offset; 
  int shiftIndex; 
  int historyIndex; 

  dfloat eps; // regularization thickness for level set function

  TimeStepper::timeStepper_t* timeStepper;
  subcell_t *subcell; 

  halo_t* traceHalo;

  dfloat *q;
  dfloat *U;
  dfloat *phi; 
  dfloat *phiH;  
  dfloat *gradq;
  dfloat *sgnq; // Depreceated
  dfloat *ssgnq; //Depreceated
  dfloat *sq; 

  dfloat *rtime; 
  occa::memory o_rtime; 


  occa::memory o_invDegree;
  // dummy delete later
  dfloat * sface; 
  occa::memory o_sface; 

  occa::memory o_gradq;
  occa::memory o_q;
  occa::memory o_U;
  occa::memory o_phi, o_phiH; 
  occa::memory o_sgnq; // Depreceated
  occa::memory o_ssgnq, o_sq; // Depreceated

  occa::memory o_Mq;

  occa::kernel advectionVolumeKernel;
  occa::kernel advectionSurfaceKernel;

  occa::kernel redistanceVolumeKernel;
  occa::kernel redistanceSurfaceKernel;

  occa::kernel gradientVolumeKernel;
  occa::kernel gradientSurfaceKernel;

  occa::kernel MassMatrixKernel;

  occa::kernel initialConditionKernel;
  occa::kernel setFlowFieldKernel;
  occa::kernel regularizedSignKernel;

  occa::kernel partialRedistanceVolumeKernel;  
  occa::kernel mixedRedistanceSurfaceKernel;  
  occa::kernel partialRedistanceSurfaceKernel; 
  occa::kernel reconstructInternalFaceKernel; 
  occa::kernel reconstructExternalFaceKernel; 
  occa::kernel reconstructFaceKernel; 

  occa::kernel projectKernel; 
  occa::kernel projectDGKernel; 
  occa::kernel reconstructKernel; 
  occa::kernel subcellComputeKernel; 
  occa::kernel indicatorMDAKernel; // Averaged Modal Decay i.e. skyline 
  occa::kernel indicatorMDHKernel; // Highest Modal Decay i.e. Pearsson and Perreira 
  // occa::kernel skyline1DKernel; 
  occa::kernel findNeighKernel; 

  occa::kernel reconstructENOKernel; 
  occa::kernel setAuxiliaryFieldKernel; 
  occa::kernel initialHistoryKernel;

  // occa::kernel subcellSignKernel; // This could be part of subcell 
  // occa::kernel subcellComputeKernel; 

  // occa::kernel subcellReconstructFaceKernel; // This could be part of subcell 
  // occa::kernel subcellSignKernel; // This could be part of subcell 
  // occa::kernel subcellComputeKernel; 
  // occa::kernel projectKernel; 
  // occa::kernel reconstructKernel; 



  lss_t() = delete;
  lss_t(mesh_t& _mesh, linAlg_t& _linAlg, settings_t& _settings):
    solver_t(_mesh, _linAlg, _settings) {}

  ~lss_t();

  //setup
  static lss_t& Setup(mesh_t& mesh, linAlg_t& linAlg,
                            lssSettings_t& settings);

  void Run();
  void Error(dfloat time, int tstep);

  void Advection(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T); 
  void Redistance(occa::memory& o_Q, occa::memory& o_RHS, const dfloat T); 
  void Redistance_Subcell(occa::memory& o_Q, occa::memory& o_sQ, 
                          occa::memory& o_RHS, occa::memory& o_sRHS,const dfloat T); 

  void Report(dfloat time, int tstep);

  void PlotFields(dfloat* Q, char *fileName);

  void rhsf(occa::memory& o_q, occa::memory& o_RHS, const dfloat dt);
  void rhsa(occa::memory& o_q, const dfloat time, const dfloat dt);

  void rhsf_subcell(occa::memory& o_Q, occa::memory &sQ, occa::memory& o_RHS,occa::memory& o_sRHS, const dfloat T);
  void reconstruct_subcell(occa::memory& o_Q, occa::memory& sQ);

  void SetupStabilizer(); 

  void DetectTroubledCells(occa::memory& o_Q, occa::memory& o_list); 
  
};


#endif

