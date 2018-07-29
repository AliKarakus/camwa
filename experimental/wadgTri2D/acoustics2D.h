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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

void acousticsSetup2D(mesh2D *mesh);

void acousticsPmlSetup2D(mesh2D *mesh,
			 dfloat xmin, dfloat xmax, // bounding box for non-pml sub-domain
			 dfloat ymin, dfloat ymax,
			 dfloat xsigma, dfloat ysigma);

void acousticsSplitPmlSetup2D(mesh2D *mesh);

void acousticsPml2D(mesh2D *mesh);

void acousticsRun2D(mesh2D *mesh);

void acousticsOccaRun2D(mesh2D *mesh);

void acousticsSplitPmlOccaRun2D(mesh2D *mesh);

void acousticsVolume2D(mesh2D *mesh);
void acousticsSurface2D(mesh2D *mesh, dfloat t);
void acousticsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void acousticsPmlUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void acousticsError2D(mesh2D *mesh, dfloat time);

void acousticsComputeVorticity2D(mesh2D *mesh, dfloat *q, int outfld, int Nfields);

void acousticsCavitySolution2D(dfloat x, dfloat y, dfloat t, dfloat *u, dfloat *v, dfloat *p);

void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t, dfloat *u, dfloat *v, dfloat *p);
