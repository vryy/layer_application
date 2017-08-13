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
lh.Print()

######################################
for str_layer in layers.layer_list:
    lh.CollapseLayer(str_layer, 1.0e-6)
print("local collapsing completed")
lh.RenumberAll()
lh.Print()
 
######################################
lh.AddGroup("group 1", ["l1", "l2"])
lh.AddGroup("group 2", ["l2", "l3"])
lh.AddGroup("group 3", ["l3", "l4"])
lh.CollapseGroup("group 1", 1.0e-6)
lh.CollapseGroup("group 2", 1.0e-6)
lh.CollapseGroup("group 3", 1.0e-6)
lh.RenumberAll()
print("global collapsing completed")
lh.Print()

lh.WriteMDPA("system")
print("WriteMDPA completed")

