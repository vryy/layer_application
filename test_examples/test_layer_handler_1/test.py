import sys
import os
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.LayerApplication import *
kernel = Kernel()

import system_layers
layers = system_layers.LayerProvider()

lh = LayerHandler()
lh.SetSpacing(0.1, 0.1, 0.1)

for str_layer in layers.layer_list:
    lh.AddLayer(str_layer, layers.layer_nodes_sets[str_layer], layers.layer_entities_sets[str_layer], layers.layer_entity_info_sets[str_layer])
gravity = Vector(3)
gravity[0] = 0.0
gravity[1] = 0.0
gravity[2] = -9.81
permeability = Matrix(3, 3)
lh["sphere"].set("DENSITY", 0.1)
lh["sphere"].set("YOUNG_MODULUS", 2000.0)
lh["sphere"].set("POISSON_RATIO", 0.3)
lh["sphere"].set("LAYER_ENTITY_TYPE", "element")
lh["sphere"].set("LAYER_ENTITY_NAME", "KinematicLinearGeo3dBezier")
lh["sphere"].set("GRAVITY", gravity)
# lh["sphere"].set("PERMEABILITY_TENSOR", permeability)

lh["surface"].set("LAYER_ENTITY_TYPE", "condition")
lh["surface"].set("LAYER_ENTITY_NAME", "FacePressureBezier")

######################################
# lh.Print()
lh.CollapseLayer("sphere", 1.0e-6)
lh.CollapseLayer("surface", 1.0e-6)
print("local collapsing completed")
lh.RenumberAll()
# lh.Print()
# 
# ######################################
lh.AddGroup("group 1", ["sphere", "surface"])
lh.CollapseGroup("group 1", 1.0e-6)
lh.RenumberAll()
print("global collapsing completed")
# lh.Print()

lh.WriteMDPA("system")
print("WriteMDPA completed")

