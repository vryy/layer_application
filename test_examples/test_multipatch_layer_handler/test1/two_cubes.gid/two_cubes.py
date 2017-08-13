##################################################################
##### KRATOS Simulation script                               #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Janosch Stascheit for TUNCONSTRUCT        #####
#####          and Giang H. Bui for SD-RUB                   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
## ATTENTION: here the order is important                    #####
##################################################################
## including kratos path                                     #####
## ATTENTION: the following lines have to be adapted to      #####
##            match your acrtual configuration               #####
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##setting up paths
kratos_libs_path = kratos_root_path+'libs' ##kratos_root/libs
kratos_applications_path = kratos_root_path+'applications' ##kratos_root/applications
##################################################################
##################################################################
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

##################################################################
##################################################################
sys.path.append('./two_cubes.gid')
parallel_type = 'shared'
if(parallel_type == 'shared'):
    import two_cubes_structural_include
    from two_cubes_structural_include import *
    model = two_cubes_structural_include.Model('two_cubes',os.getcwd()+"/")
elif(parallel_type == 'distributed'):
    import two_cubes_structural_parallel_include
    from two_cubes_structural_parallel_include import *
    model = two_cubes_structural_parallel_include.Model('two_cubes',os.getcwd()+"/")
model.InitializeModel()

##################################################################
###  SIMULATION  #################################################
##################################################################

