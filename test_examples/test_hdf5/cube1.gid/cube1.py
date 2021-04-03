##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
#####  copyright (c) (2009, 2010, 2011, 2012, 2013)          #####
#####   by CIMNE, Barcelona, Spain and Janosch Stascheit     #####
#####           for TUNCONSTRUCT                             #####
#####  and (c) 2014, 2015, 2016, 2017, 2018, 2019            #####
#####     by Hoang-Giang Bui for SFB837                      #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## This file is generated on Sa 14. Mar 00:15:32 CET 2020 
##################################################################
import sys
import os
import math
import time as time_module
##################################################################
##################################################################
sys.path.append('./cube1.gid')
import cube1_include
from cube1_include import *
# calculate insitu-stress for geology_virgin.gid
model = cube1_include.Model('cube1',os.getcwd()+"/",os.getcwd()+"/")
model.InitializeModel()
##################################################################
###  SIMULATION  #################################################
start_time = time_module.time()
##################################################################

tol = 1.0e-6
prescribed_nodes = []
for node in model.model_part.Nodes:
    if abs(node.X0) < tol:
        node.Fix(DISPLACEMENT_X)
    if abs(node.Y0) < tol:
        node.Fix(DISPLACEMENT_Y)
    if abs(node.Z0) < tol:
        node.Fix(DISPLACEMENT_Z)
    if abs(node.Z0 - 1.0) < tol:
        node.Fix(DISPLACEMENT_Z)
        prescribed_nodes.append(node)

time = 0.0
model.Solve(time, 0, 0, 0, 0)

du = -0.1
for node in prescribed_nodes:
    node.SetSolutionStepValue(DISPLACEMENT_Z, du)

time = 1.0
model.Solve(time, 0, 0, 0, 0)
model.WriteOutput(time)

for node in model.model_part.Nodes:
    print("node " + str(node.Id) + " DISPLACEMENT: " + str(node.GetSolutionStepValue(DISPLACEMENT)))

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
