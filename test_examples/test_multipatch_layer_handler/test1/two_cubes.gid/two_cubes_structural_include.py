##################################################################
##################### shared_include.py   ########################
##################################################################
##### KRATOS Multiphysics                                    #####
##### include file for shared-memory simulation              #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MKLSolversApplication import *
from KratosMultiphysics.LayerApplication import *
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path ):
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("Kratos_simulation")
        self.path = path
        self.problem_name = problem_name
        ##################################################################
        ## ADD VARIABLES #################################################
        ##################################################################
        import structural_solver_static
        structural_solver_static.AddVariables(self.model_part)

        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        # defining linear solver
        self.structure_linear_solver = MKLPardisoSolver()

        # defining builder_and_solver
        self.builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)

        # defining time scheme
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # defining convergence criteria
        self.conv_criteria = MultiPhaseFlowCriteria(       1e-06,        1e-09)

        # defining solving strategy flags
        self.MaxNewtonRapshonIterations = 30
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True

        # defining solving strategy
        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part, \
                                                         self.time_scheme, \
                                                         self.structure_linear_solver, \
                                                         self.conv_criteria, \
                                                         self.builder_and_solver, \
                                                         self.MaxNewtonRapshonIterations, \
                                                         self.CalculateReactionFlag, \
                                                         self.ReformDofSetAtEachStep, \
                                                         self.MoveMeshFlag)

        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.model_part_io.ReadModelPart(self.model_part)
        self.meshWritten = False
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        self.element_assignments = {}
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
                self.element_assignments[int(val_set[1])] = int(val_set[2])
        self.cond_file.close()

        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)


        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
        structural_solver_static.AddDofs(self.model_part)

    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False

        self.model_part.Properties[1].SetValue(DENSITY,         7620 )
        self.model_part.Properties[1].SetValue(YOUNG_MODULUS,      2.1e+05 )
        self.model_part.Properties[1].SetValue(POISSON_RATIO,          0.3 )
        self.model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
        print("Linear elastic (Saint-Venant Kirschoff) is set")

        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        model_layers = __import__(self.problem_name+"_layers")
        ## ELEMENTS on layers ############################################
        self.layer_sets = model_layers.ReadLayerSets()
        ## NODES on layers ###############################################
        self.layer_nodes_sets = model_layers.ReadLayerNodesSets()
        ##################################################################
        print "layer sets stored"
        ##################################################################
        ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        ##################################################################
        self.node_groups = model_layers.ReadNodeGroups()
        print "node groups stored"
        ##################################################################
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
        print "activation utility initialized"
        self.model_part.Check( self.model_part.ProcessInfo )
        print "model successfully initialized"

    def WriteOutput( self, time ):
        self.gid_io.InitializeMesh( time )
        mesh = self.model_part.GetMesh()
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        print("mesh written...")
        self.gid_io.FinalizeMesh()
        self.gid_io.InitializeResults( time, self.model_part.GetMesh() )
        print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.Nodes, time, 0)
        self.gid_io.FinalizeResults()

    def WriteRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time)
        serializer = Serializer(fn)
        serializer.Save("ModelPart", self.model_part)
        serializer = 0
        print("Write restart data to " + fn + ".rest completed")

    def LoadRestartFile( self, time ):
        fn = self.problem_name + "_" + str(time)
        serializer = Serializer(fn)
        serializer.Load("ModelPart", self.model_part)
        serializer = 0
        print("Load restart data from " + fn + ".rest completed")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
##################################################################
