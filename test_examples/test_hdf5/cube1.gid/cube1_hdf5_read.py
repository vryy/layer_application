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

time = 1.0
filename = model.problem_name + "_" + str(time) + ".h5"
hdf5_post = HDF5PostUtility(model.problem_name + "_" + str(time) + ".h5", "Read-Only")
hdf5_post.ReadNodalResults(DISPLACEMENT, model.model_part)
print("Read hdf5 completed")

for node in model.model_part.Nodes:
    print("node " + str(node.Id) + " DISPLACEMENT: " + str(node.GetSolutionStepValue(DISPLACEMENT)))

## Reference results:
# node 1 DISPLACEMENT: [3](0,0.03,-0.1)
# node 2 DISPLACEMENT: [3](0,0.03,0)
# node 3 DISPLACEMENT: [3](0,0,-0.1)
# node 4 DISPLACEMENT: [3](0.03,0.03,-0.1)
# node 5 DISPLACEMENT: [3](0,0,0)
# node 6 DISPLACEMENT: [3](0.03,0.03,0)
# node 7 DISPLACEMENT: [3](0.03,0,-0.1)
# node 8 DISPLACEMENT: [3](0.03,0,0)

##################################################################
###  END OF SIMULATION  ##########################################
end_time = time_module.time()
print("Calculation time: " + str(end_time - start_time) + " s")
timer = Timer()
print(timer)
##################################################################
