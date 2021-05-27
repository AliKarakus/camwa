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
#define CIRCLE_TEST 0
#define SQUARE_TEST 1
#define INTERSECTINGCIRCLE_TEST 0


#include "lss.hpp"

void lss_t::Error(dfloat time, int tstep){

   dlong   Nentries = mesh.Nelements*mesh.Np;
   dfloat *test_0   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
   dfloat *test_1   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
   occa::memory o_test = mesh.device.malloc(Nentries*sizeof(dfloat), test_0);

 
   o_test.copyFrom(o_phi); 
  
#if 1
   mesh.ogs->GatherScatter(o_test, ogs_dfloat, ogs_add, ogs_sym);
   linAlg.amx(Nentries, 1.0, o_invDegree, o_test); 
#endif
   o_test.copyTo(phi);


  for(int n=0; n<Nentries; n++){
    const dfloat xn = mesh.x[n]; 
    const dfloat yn = mesh.y[n]; 

    #if CIRCLE_TEST
    const dfloat pn    = phi[n];
    const dfloat exact = sqrt(xn*xn + yn*yn) - 1.0;
    
    test_0[n] = fabs(pn - exact); 
    test_1[n] = 0.0; 
    
    const dfloat beps = 0.25; 
    if(fabs(exact)<= beps )
      test_1[n] = fabs(pn-exact); 
   
    #elif SQUARE_TEST
     const dfloat dx = fabs(xn) - 1.0; 
     const dfloat dy = fabs(yn) - 1.0; 
     //
     const dfloat tdx = std::max(dx,0.0); 
     const dfloat tdy = std::max(dy,0.0); 
     const dfloat d  = sqrt(tdx*tdx + tdy*tdy); 

     const dfloat pn    = phi[n];
     const dfloat exact = d + std::min(std::max(dx,dy),0.0);
    
     test_0[n] = fabs(pn - exact); 
     test_1[n] = 0.0; 
    
    const dfloat beps = 0.25; 
    if(fabs(exact)<= beps )
      test_1[n] = fabs(pn-exact); 

    #elif INTERSECTINGCIRCLE_TEST



    #endif 
  }



  o_test.copyFrom(test_0); 
  MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
  dfloat normg2 = sqrt(linAlg.innerProd(Nentries, o_test, o_Mq, comm));

  o_test.copyFrom(test_1); 
  MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
  dfloat norml2 = sqrt(linAlg.innerProd(Nentries, o_test, o_Mq, comm));


   if(mesh.rank==0){
    printf("%5.2f (%d), %.8e %.8e (time, timestep, norm_g, norm_l)\n", time, tstep, normg2, norml2);    
   }

   free(test_0); free(test_1);


  // if(mesh.rank==0)
  //   printf("%5.2f (%d), %.8e (time, timestep, norm)\n", time, tstep, norm2);

  // if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {

  //   o_q.copyTo(q);
  //   // o_phi.copyTo(phi); 
  //   if(redistance){
  //     // o_sgnq.copyTo(sgnq);
  //     subcell->o_ElementList.copyTo(subcell->ElementList);
  //   }   
  //   // output field files
  //   string name;
  //   settings.getSetting("OUTPUT FILE NAME", name);
  //   char fname[BUFSIZ];
  //   sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);

  //   PlotFields(q, fname);
  // }
}
