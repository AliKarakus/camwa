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
#define CIRCLE_TEST 1
#define SQUARE_TEST 0
#define INTERSECTINGCIRCLE_TEST 0


#include "lss.hpp"

void lss_t::Report(dfloat time, int tstep){

  // Error(time, tstep);

  static int frame=0;

  if (settings.compareSetting("OUTPUT TO FILE","TRUE")) {
    dlong   Nentries = mesh.Nelements*mesh.Np;
    dfloat *test   = (dfloat *)calloc(Nentries, sizeof(dfloat)); 
    occa::memory o_test = mesh.device.malloc(Nentries*sizeof(dfloat), test);
    
    o_test.copyFrom(o_phi);


    MassMatrixKernel(mesh.Nelements, mesh.o_ggeo, mesh.o_MM, o_test, o_Mq);
    dfloat norml2 = sqrt(linAlg.innerProd(Nentries, o_test, o_Mq, comm));

    if(mesh.rank==0){
        printf("%5.2f, %.4e (time,  L2_norm)\n", 
          time,  norml2);  
      }


    mesh.ogs->GatherScatter(o_test, ogs_dfloat, ogs_add, ogs_sym);
    linAlg.amx(Nentries, 1.0, o_invDegree, o_test); 
    o_test.copyTo(phi);

    if(redistance){
      subcell->o_ElementList.copyTo(subcell->ElementList);
    }   
    // output field files
    string name;
    settings.getSetting("OUTPUT FILE NAME", name);
        
    if(frame==0){
      char fname_tri[BUFSIZ];
      sprintf(fname_tri, "%s_tri.dat", name.c_str());
      writeConnectivity(q,fname_tri);
    }

    char fname_data[BUFSIZ];
    sprintf(fname_data, "%s_%d.dat", name.c_str(), frame++);
    writeData(q, fname_data);
  }


  //   o_q.copyTo(q);

  //   if(redistance){
  //     subcell->o_ElementList.copyTo(subcell->ElementList);
  //   }   
  //   // output field files
  //   string name;
  //   settings.getSetting("OUTPUT FILE NAME", name);
  //   char fname[BUFSIZ];
  //   // sprintf(fname, "%s_%04d_%04d.vtu", name.c_str(), mesh.rank, frame++);
  //   sprintf(fname, "%s_%d.vtu", name.c_str(), frame++);

  //   PlotFields(q, fname);
  // }
}
