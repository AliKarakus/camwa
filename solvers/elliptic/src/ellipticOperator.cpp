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

#include "elliptic.hpp"

void elliptic_t::Operator(occa::memory &o_q, occa::memory &o_Aq){

  if(continuous){
    int mapType = (mesh.elementType==HEXAHEDRA &&
                   options.compareArgs("ELEMENT MAP", "TRILINEAR")) ? 1:0;

    // int integrationType = (mesh.elementType==HEXAHEDRA &&
    //                        options.compareArgs("ELLIPTIC INTEGRATION", "CUBATURE")) ? 1:0;

    if(mesh.NglobalGatherElements) {

      // if(integrationType==0) { // GLL or non-hex
        if(mapType==0)
          partialAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList,
                          mesh.o_ggeo, mesh.o_Dmatrices, mesh.o_Smatrices, mesh.o_MM, lambda, o_q, o_Aq);
        else
          partialAxKernel(mesh.NglobalGatherElements, mesh.o_globalGatherElementList,
                          elliptic->o_EXYZ, elliptic->o_gllzw, mesh.o_Dmatrices, mesh.o_Smatrices, mesh.o_MM, lambda, o_q, o_Aq);
      // } else {
      //   partialCubatureAxKernel(mesh.NglobalGatherElements,
      //                           mesh.o_globalGatherElementList,
      //                           mesh.o_cubggeo,
      //                           mesh.o_cubD,
      //                           mesh.o_cubInterpT,
      //                           lambda, o_q, o_Aq);
      // }
    }

    ogsGatherScatterStart(o_Aq, ogsDfloat, ogsAdd, ogs);

    if(mesh.NlocalGatherElements){
      // if(integrationType==0) { // GLL or non-hex
        if(mapType==0)
          partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
                          mesh.o_ggeo, mesh.o_Dmatrices, mesh.o_Smatrices, mesh.o_MM, lambda, o_q, o_Aq);
        else
          partialAxKernel(mesh.NlocalGatherElements, mesh.o_localGatherElementList,
                          elliptic->o_EXYZ, elliptic->o_gllzw, mesh.o_Dmatrices, mesh.o_Smatrices, mesh.o_MM, lambda, o_q, o_Aq);
      // } else {
      //   partialCubatureAxKernel(mesh.NlocalGatherElements,
      //                           mesh.o_localGatherElementList,
      //                           mesh.o_cubggeo,
      //                           mesh.o_cubD,
      //                           mesh.o_cubInterpT,
      //                           lambda,
      //                           o_q,
      //                           o_Aq);
      // }
    }

    // finalize gather using local and global contributions
    ogsGatherScatterFinish(o_Aq, ogsDfloat, ogsAdd, ogs);

#if USE_NULL_BOOST==1
    if(mesh.allNeumann) {
      dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
      dfloat alpha = linAlg->innerProd(Nentries, o_invDegree, o_q, comm);
      alpha *= allNeumannPenalty*allNeumannScale*allNeumannScale;

      linAlg->add(Nentries, alpha, o_Aq);
    }
#endif

    //post-mask
    if (Nmasked)
      maskKernel(Nmasked, o_maskIds, o_Aq);

  } else if(ipdg) {

    int Nentries = mesh.Np;
    mesh.HaloExchangeStart(o_q, Nentries, ogsDfloat);

    if(mesh.Nelements) {
      dlong offset = 0;
      partialGradientKernel(mesh.Nelements,
                            offset,
                            mesh.o_vgeo,
                            mesh.o_Dmatrices,
                            o_q,
                            o_grad);
    }

    if(mesh.NinternalElements)
      partialIpdgKernel(mesh.NinternalElements,
                        mesh.o_internalElementIds,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        lambda,
                        tau,
                        mesh.o_vgeo,
                        mesh.o_sgeo,
                        o_EToB,
                        mesh.o_Dmatrices,
                        mesh.o_LIFTT,
                        mesh.o_MM,
                        o_grad,
                        o_Aq);

    mesh.HaloExchangeFinish(o_q, Nentries, ogsDfloat);

    if(mesh.totalHaloPairs){
      dlong offset = mesh.Nelements;
      partialGradientKernel(mesh.totalHaloPairs,
                            offset,
                            mesh.o_vgeo,
                            mesh.o_Dmatrices,
                            o_q,
                            o_grad);
    }

    if(mesh.NnotInternalElements) {
      partialIpdgKernel(mesh.NnotInternalElements,
                        mesh.o_notInternalElementIds,
                        mesh.o_vmapM,
                        mesh.o_vmapP,
                        lambda,
                        tau,
                        mesh.o_vgeo,
                        mesh.o_sgeo,
                        o_EToB,
                        mesh.o_Dmatrices,
                        mesh.o_LIFTT,
                        mesh.o_MM,
                        o_grad,
                        o_Aq);
    }

#if USE_NULL_BOOST==1
    //rank 1 augmentation if all BCs are Neumann
    //TODO this could probably be moved inside the Ax kernel for better performance
    if(allNeumann) {
      dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
      alpha = linAlg->sum(Nentries, o_q, comm);
      alpha *= allNeumannPenalty*allNeumannScale*allNeumannScale;

      linAlg->add(Nentries, alpha, o_Aq);
    }
#endif
  }
}
