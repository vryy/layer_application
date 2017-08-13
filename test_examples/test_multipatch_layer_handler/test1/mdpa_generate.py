import sys
import os
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.LayerApplication import *
kernel = Kernel()

####################################################################
#### mdpa generation for system model
####################################################################

## import layer data
import system_layers
layers = system_layers.LayerProvider()

## define the layer handler
lh = MultipatchLayerHandler()
for str_layer in layers.layer_list:
    lh.AddLayer(str_layer, layers.layer_nodes_sets[str_layer], layers.layer_entities_sets[str_layer], layers.layer_entity_info_sets[str_layer])

## add the patch boundary
lh['upper'].AddTable('bottom', layers.layer_boundary_marker['upper']['bottom'])
lh['lower'].AddTable('top', layers.layer_boundary_marker['lower']['top'])

## add the patch boundary
lh.AddLayerConnection('upper', 'bottom', 'lower', 'top')

## define the attributes of each layer
lh['upper'].set("LAYER_ENTITY_TYPE", "element")
lh['upper'].set("LAYER_ENTITY_NAME", "KinematicLinear3D8N")

lh['lower'].set("LAYER_ENTITY_TYPE", "element")
lh['lower'].set("LAYER_ENTITY_NAME", "KinematicLinear3D8N")

## renumber all the nodes of the system, with account to patches
lh.Rebuild()
#lh.RenumberAll()

## here we write the nodes and entities to the mdpa
lh.WriteMDPA("twocubes")
print("WriteMDPA completed for system")

## export the layer file
lh.WriteLayers("twocubes_layers")
print("WriteLayers completed for system")

