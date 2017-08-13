# geometry layers data file for isogeometric application
# (c) 2016 Hoang Giang Bui, Ruhr-University Bochum
# This file is created at 05-May-2016 09:41:57

class LayerProvider():
    def __init__(self):
        self.layer_list = []
        self.layer_attributes = {}
        self.layer_nodes_sets = {}
        self.layer_entities_sets = {}
        self.layer_entity_info_sets = {}
        self.layer_boundary_marker = {}

        ## begin layer_info for layer upper
        self.layer_list.append('upper')

        # nodal info
        current_nodal_set = {}
        current_nodal_set[1] = [   0.00000000,  1.00000000,  2.00000000]
        current_nodal_set[2] = [   0.00000000,  1.00000000,  1.00000000]
        current_nodal_set[3] = [   0.00000000,  0.00000000,  2.00000000]
        current_nodal_set[4] = [   1.00000000,  1.00000000,  2.00000000]
        current_nodal_set[5] = [   0.00000000,  0.00000000,  1.00000000]
        current_nodal_set[6] = [   1.00000000,  1.00000000, 1.00000000]
        current_nodal_set[7] = [   1.00000000,  0.00000000,  2.00000000]
        current_nodal_set[8] = [   1.00000000,  0.00000000,  1.00000000]

        self.layer_nodes_sets['upper'] = current_nodal_set

        # entity connectivities
        current_entity_set = {}
        current_entity_set[2] = [7, 4, 6, 8, 3, 1, 2, 5]
        self.layer_entities_sets['upper'] = current_entity_set

        # entity data
        current_entity_info = {}
        self.layer_entity_info_sets['upper'] = current_entity_info

        # layer boundary marker
        boundary_marker = {}
        boundary_marker['bottom'] = [5, 8, 6, 2]
        self.layer_boundary_marker['upper'] = boundary_marker

        # layer attributes
        self.layer_attributes['ground'] = {}
        ################### end layer_info upper ###################

        ## begin layer_info for layer lower
        self.layer_list.append('lower')

        # nodal info
        current_nodal_set = {}
        current_nodal_set[1] = [   0.00000000,  1.00000000,  1.00000000]
        current_nodal_set[2] = [   0.00000000,  0.00000000,  1.00000000]
        current_nodal_set[3] = [   1.00000000,  1.00000000,  1.00000000]
        current_nodal_set[4] = [   1.00000000,  0.00000000,  1.00000000]
        current_nodal_set[5] = [   0.00000000,  1.00000000,  0.00000000]
        current_nodal_set[6] = [   0.00000000,  0.00000000,  0.00000000]
        current_nodal_set[7] = [   1.00000000,  1.00000000,  0.00000000]
        current_nodal_set[8] = [   1.00000000,  0.00000000,  0.00000000]
        self.layer_nodes_sets['lower'] = current_nodal_set

        # entity connectivities
        current_entity_set = {}
        current_entity_set[1] = [4, 3, 7, 8, 2, 1, 5, 6]
        self.layer_entities_sets['lower'] = current_entity_set

        # entity data
        current_entity_info = {}
        self.layer_entity_info_sets['lower'] = current_entity_info

        # layer boundary marker
        boundary_marker = {}
        boundary_marker['top'] = [2, 4, 3, 1]
        self.layer_boundary_marker['lower'] = boundary_marker

        # layer attributes
        self.layer_attributes['lower'] = {}
        ################### end layer_info excavation ###################

