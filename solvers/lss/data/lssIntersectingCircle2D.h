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


// Level-Set function
#define lssInitialConditions2D(t, x, y, q) \
{                                       \
 const dfloat rc = 1.0;          \
  const dfloat ac = 0.7;          \
  const dfloat test1 = (ac-x)/sqrt((ac-x)*(ac-x) + y*y); \
  const dfloat test2 = (ac+x)/sqrt((ac+x)*(ac+x) + y*y); \
  const dfloat scale = pow((x-rc) * (x-rc) + (y-rc)*(y-rc)+0.1,1.0); \
if( ( test1 >= ac/rc ) && ( test2 >= ac/rc ) ){ \
*(q) = - scale* min(sqrt( x*x+  ( y-sqrt( rc*rc-ac*ac ) )*( y-sqrt( rc*rc-ac*ac))), sqrt( x*x+( y+sqrt( rc*rc-ac*ac))*(y+sqrt( rc*rc-ac*ac)))); \
} \
else{ \
*(q)= scale*min( sqrt( (x+ac)*(x+ac)+ y*y)-rc, sqrt( (x-ac)*(x-ac)+ y*y)-rc );\
} \
}

// Level-Set function
#define lssExactSolution2D(t, x, y, q) \
{                                       \
 const dfloat rc = 1.0;          \
  const dfloat ac = 0.7;          \
  const dfloat test1 = (ac-x)/sqrt((ac-x)*(ac-x) + y*y); \
  const dfloat test2 = (ac+x)/sqrt((ac+x)*(ac+x) + y*y); \
  const dfloat scale = 1.0 \
if( ( test1 >= ac/rc ) && ( test2 >= ac/rc ) ){ \
*(q) = - scale* min(sqrt( x*x+  ( y-sqrt( rc*rc-ac*ac ) )*( y-sqrt( rc*rc-ac*ac))), sqrt( x*x+( y+sqrt( rc*rc-ac*ac))*(y+sqrt( rc*rc-ac*ac)))); \
} \
else{ \
*(q)= scale*min( sqrt( (x+ac)*(x+ac)+ y*y)-rc, sqrt( (x-ac)*(x-ac)+ y*y)-rc );\
} \
}

 
// LS Advective field
#define lssAdvectionField2D(t, x, y, q, u, v) \
{                                       \
  const dfloat PERIOD = 8.f; 			\
  const dfloat xs = x + 0.5; 			\
  const dfloat ys = y + 0.5; 			\
 *(u) = sin(M_PI*xs)*sin(M_PI*xs)*sin(2.f*M_PI*ys)*cos(M_PI*t/PERIOD); \
 *(v) =-sin(M_PI*ys)*sin(M_PI*ys)*sin(2.f*M_PI*xs)*cos(M_PI*t/PERIOD); \
}


// Boundary conditions
/* wall 1, outflow 2 */
#define lssDirichletConditions2D(bc, t, x, y, nx, ny, qM, qB) \
{                                       \
  if(bc==1){                            \
    *(qB) = qM;                        \
  } else if(bc==2){                     \
    *(qB) = qM;                         \
  }                                     \
}