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

import structural_solver_advanced

time = 1
ref_file_name = model.problem_name+'_'+str(time)+'.post.bin'
# print(ref_file_name)
ref_model_part = ModelPart(ref_file_name)
structural_solver_advanced.AddVariables(ref_model_part)
reader = GiDBinaryReader(ref_file_name)
print(reader)
# mesh_list = reader.GetMeshesName()
# print("mesh_list:", mesh_list)
disp_list = reader.ReadNodalArray1DValues(DISPLACEMENT, 0)
print("disp_list:", disp_list)
temp_list = reader.ReadNodalScalarValues(TEMPERATURE, 0)
print("temp_list:", temp_list)
# reader.ReadMesh(ref_model_part)
# reader.ReadNodalDisplacements(ref_model_part)
# print(ref_model_part)
# del reader

# for node in ref_model_part.Nodes:
#     print("node " + str(node.Id) + " DISPLACEMENT: " + str(node.GetSolutionStepValue(DISPLACEMENT)))

## reference results:
# node 1 DISPLACEMENT: [3](0,0.03,-0.1)
# node 2 DISPLACEMENT: [3](0,0.03,0)
# node 3 DISPLACEMENT: [3](0,0,-0.1)
# node 4 DISPLACEMENT: [3](0.03,0.03,-0.1)
# node 5 DISPLACEMENT: [3](0,0,0)
# node 6 DISPLACEMENT: [3](0.03,0.03,0)
# node 7 DISPLACEMENT: [3](0.03,0,-0.1)
# node 8 DISPLACEMENT: [3](0.03,0,0)

# node 1 TEMPERATURE: 1.1
# node 2 TEMPERATURE: 2.1
# node 3 TEMPERATURE: 3.1
# node 4 TEMPERATURE: 4.1
# node 5 TEMPERATURE: 5.1
# node 6 TEMPERATURE: 6.1
# node 7 TEMPERATURE: 7.1
# node 8 TEMPERATURE: 8.1

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
