#!/usr/bin/env python3

#####################################################################################
#
#The MIT License (MIT)
#
#Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#####################################################################################

from test import *

lssCircle2D = lssDir + "/data/lssCircle2D.h"
lssElliptic2D = lssDir + "/data/lssElliptic2D.h"
lssIntersecting2D = lssDir + "/data/lssIntersectingCircle2D.h"
lssSquare2D = lssDir + "/data/lssSquare2D.h"
lssMultiple2D = lssDir + "/data/lssMultipleInterfaces2D.h"

def lssSettings(rcformat="2.0", data_file=lssCircle2D, redist="TRUE", advect="FALSE", advect_tpe="COLLOCATION",
               mesh="circle_2.msh", dim=2, element=3, nx=10, ny=10, nz=10, 
               boundary_flag=1, box_dimx=4, box_dimy=4,
               degree=5, stab_type="SUBCELL", subcell_no=6, indicator_type="MDA", minor_grid="WARPBLEND",time_reconst="ENO2"
               ,thread_model="CUDA", platform_number=0, device_number=0,  
               time_integrator="LSERK4", cfl=1.0, start_time=0.0, final_time=2.5,
               output_to_file="TRUE", out_interval=10, out_file_name="test1"):
  return [setting_t("FORMAT", rcformat),
          setting_t("DATA FILE", data_file),
          setting_t("MESH FILE", mesh),
          setting_t("MESH DIMENSION", dim),
          setting_t("ELEMENT TYPE", element),
          setting_t("BOX NX", nx),
          setting_t("BOX NY", ny),
          setting_t("BOX NZ", nz),
          setting_t("BOX BOUNDARY FLAG", boundary_flag),
          setting_t("POLYNOMIAL DEGREE", degree),
          setting_t("THREAD MODEL", thread_model),
          setting_t("PLATFORM NUMBER", platform_number),
          setting_t("ADVECTION SOLVER", advect),
          setting_t("REDISTANCE SOLVER", redist),
          setting_t("ADVECTION TYPE", advect_tpe),
          setting_t("STABILIZATION", stab_type),
          setting_t("SUBCELL NUMBER", subcell_no),
          setting_t("INDICATOR TYPE", indicator_type),
          setting_t("TIME RECONSTRUCTION", time_reconst),
          setting_t("SUBCELL MINOR GRID", minor_grid),
          setting_t("TIME INTEGRATOR", time_integrator),
          # setting_t("CFL NUMBER", cfl),
          setting_t("START TIME", start_time),
          setting_t("FINAL TIME", final_time),
          setting_t("OUTPUT TO FILE", output_to_file),
          setting_t("OUTPUT FILE NAME", out_file_name),
          setting_t("OUTPUT INTERVAL", out_interval)]

def main():
  failCount=0;

  failCount += test(name="testSmoothCircle_h1",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssCircle2D,dim=2,final_time=2.0, mesh="circle_1.msh", degree=3,subcell_no=4,
                                         out_file_name="circle_h1"),
                    referenceNorm=3.11475459668639)

  createContour(triFile="circle_h1_tri.dat", dataFile="circle_h1_1.dat", outFile="smoothCircle_b.pdf", levelmin=-0.9, levelmax=0.9,leveli=0.1)


  failCount += test(name="testSmoothCircle_h2",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssCircle2D,dim=2,final_time=2.0, mesh="circle_2.msh", degree=3,subcell_no=4,
                                         out_file_name="circle_h2"),
                    referenceNorm=3.11351985632056)

  createContour(triFile="circle_h2_tri.dat", dataFile="circle_h2_0.dat", outFile="smoothCircle_a.pdf", levelmin=-0.9, levelmax=0.9,leveli=0.1)
  createContour(triFile="circle_h2_tri.dat", dataFile="circle_h2_1.dat", outFile="smoothCircle_c.pdf", levelmin=-0.9, levelmax=0.9,leveli=0.1)

  failCount += test(name="testElliptic_h2",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssElliptic2D,dim=2,final_time=2.5,mesh="circle_2.msh",
                                         out_file_name="elliptic_h2"),
                    referenceNorm=3.95462790564238)

  createContour(triFile="elliptic_h2_tri.dat", dataFile="elliptic_h2_1.dat", outFile="elliptic_b.pdf", levelmin=-0.4, levelmax=1.5,leveli=0.1)

  failCount += test(name="testElliptic_h4",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssElliptic2D,final_time=2.5, dim=2,mesh="circle_4.msh",
                                         out_file_name="elliptic_h4"),
                    referenceNorm=3.95461727175103)

  createContour(triFile="elliptic_h4_tri.dat", dataFile="elliptic_h4_0.dat", outFile="elliptic_a.pdf", levelmin=-0.4, levelmax=1.5,leveli=0.1)
  createContour(triFile="elliptic_h4_tri.dat", dataFile="elliptic_h4_1.dat", outFile="elliptic_c.pdf", levelmin=-0.4, levelmax=1.5,leveli=0.1)


  failCount += test(name="testIntersecting_h1",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssIntersecting2D,dim=2,final_time=1.5, mesh="circle_1.msh",
                                         out_file_name="intersectingCircles_h1"),
                    referenceNorm=2.27903550574234)

  createContour(triFile="intersectingCircles_h1_tri.dat", dataFile="intersectingCircles_h1_1.dat", outFile="intersectingCircles_b.pdf", levelmin=-0.9, levelmax=0.9,leveli=0.1)

  failCount += test(name="testIntersecting_h2",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssIntersecting2D,dim=2,final_time=1.5,mesh="circle_2.msh",
                                         out_file_name="intersectingCircles_h2"),
                    referenceNorm=2.27807224946598)

  createContour(triFile="intersectingCircles_h2_tri.dat", dataFile="intersectingCircles_h2_0.dat", outFile="intersectingCircles_a.pdf", levelmin=-0.9, levelmax=0.9,leveli=0.1)
  createContour(triFile="intersectingCircles_h2_tri.dat", dataFile="intersectingCircles_h2_1.dat", outFile="intersectingCircles_c.pdf", levelmin=-0.9, levelmax=0.9,leveli=0.1)


  failCount += test(name="testSquare_h2",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssSquare2D,dim=2,final_time=1.5, mesh="circle_2.msh",
                                         out_file_name="square_h2"),
                    referenceNorm=2.45461854838549)

  createContour(triFile="square_h2_tri.dat", dataFile="square_h2_1.dat", outFile="square_b.pdf", levelmin=-0.9, levelmax=1.0,leveli=0.1)

  failCount += test(name="testSquare_h4",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssSquare2D,dim=2,final_time=1.5,mesh="circle_4.msh",
                                         out_file_name="square_h4"),
                    referenceNorm=2.45184829788086)

  createContour(triFile="square_h4_tri.dat", dataFile="square_h4_0.dat", outFile="square_a.pdf", levelmin=-0.9, levelmax=1.0,leveli=0.1)
  createContour(triFile="square_h4_tri.dat", dataFile="square_h4_1.dat", outFile="square_c.pdf", levelmin=-0.9, levelmax=1.0,leveli=0.1)

  failCount += test(name="testMultiple_h4",
                    cmd=lssBin,
                    settings=lssSettings(element=3,data_file=lssMultiple2D,dim=2,final_time=1.5, mesh="circle_4.msh",
                                         out_file_name="multiple_h4"),
                    referenceNorm=1.7077622642489)

  createContour(triFile="multiple_h4_tri.dat", dataFile="multiple_h4_0.dat", outFile="multiple_a.pdf", levelmin=-0.5, levelmax=1.0,leveli=0.05)
  createContour(triFile="multiple_h4_tri.dat", dataFile="multiple_h4_1.dat", outFile="multiple_b.pdf", levelmin=-0.5, levelmax=1.0,leveli=0.05)

  # #clean up
  # for file_name in os.listdir(testDir):
  #   if file_name.endswith('.vtu'):
  #     os.remove(testDir + "/" + file_name)

  return failCount

if __name__ == "__main__":
  failCount=0;
  failCount+=main()
  sys.exit(failCount)
