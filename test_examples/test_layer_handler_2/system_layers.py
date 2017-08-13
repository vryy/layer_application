# geometry layers data file for isogeometric application
# (c) 2014 Hoang Giang Bui, Ruhr-University Bochum

class LayerProvider():
	def __init__(self):
		self.layer_list = []
		self.layer_attributes = {}
		self.layer_nodes_sets = {}
		self.layer_entities_sets = {}
		self.layer_entity_info_sets = {}

		## begin layer_info for layer l1
		self.layer_list.append('l1')

		# nodal info
		current_nodal_set = {}
		current_nodal_set[1] = [0.000000, 0.000000, 0.000000]
		current_nodal_set[2] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[3] = [1.000000, 1.000000, 0.000000]
		current_nodal_set[4] = [0.000000, 1.000000, 0.000000]
		self.layer_nodes_sets['l1'] = current_nodal_set

		# entity connectivities
		current_entity_set = {}
		current_entity_set[1] = [1, 2, 3, 4]
		self.layer_entities_sets['l1'] = current_entity_set

		# entity data
		current_entity_info = {}
		self.layer_entity_info_sets['l1'] = current_entity_info

		# layer attributes
		self.layer_attributes['l1'] = {}
		################### end layer_info l1 ###################

		## begin layer_info for layer l2
		self.layer_list.append('l2')

		# nodal info
		current_nodal_set = {}
		current_nodal_set[1] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[2] = [2.000000, 0.000000, 0.000000]
		current_nodal_set[3] = [2.000000, 1.000000, 0.000000]
		current_nodal_set[4] = [1.000000, 1.000000, 0.000000]
		self.layer_nodes_sets['l2'] = current_nodal_set

		# entity connectivities
		current_entity_set = {}
		current_entity_set[1] = [1, 2, 3, 4]
		self.layer_entities_sets['l2'] = current_entity_set

		# entity data
		current_entity_info = {}
		self.layer_entity_info_sets['l2'] = current_entity_info

		# layer attributes
		self.layer_attributes['l2'] = {}
		################### end layer_info l2 ###################

		## begin layer_info for layer l3
		self.layer_list.append('l3')

		# nodal info
		current_nodal_set = {}
		current_nodal_set[1] = [1.000000, 1.000000, 0.000000]
		current_nodal_set[2] = [2.000000, 1.000000, 0.000000]
		current_nodal_set[3] = [2.000000, 2.000000, 0.000000]
		current_nodal_set[4] = [1.000000, 2.000000, 0.000000]
		self.layer_nodes_sets['l3'] = current_nodal_set

		# entity connectivities
		current_entity_set = {}
		current_entity_set[1] = [1, 2, 3, 4]
		self.layer_entities_sets['l3'] = current_entity_set

		# entity data
		current_entity_info = {}
		self.layer_entity_info_sets['l3'] = current_entity_info

		# layer attributes
		self.layer_attributes['l3'] = {}
		################### end layer_info l3 ###################

		## begin layer_info for layer l4
		self.layer_list.append('l4')

		# nodal info
		current_nodal_set = {}
		current_nodal_set[1] = [0.000000, 1.000000, 0.000000]
		current_nodal_set[2] = [1.000000, 1.000000, 0.000000]
		current_nodal_set[3] = [1.000000, 2.000000, 0.000000]
		current_nodal_set[4] = [0.000000, 2.000000, 0.000000]
		self.layer_nodes_sets['l4'] = current_nodal_set

		# entity connectivities
		current_entity_set = {}
		current_entity_set[1] = [1, 2, 3, 4]
		self.layer_entities_sets['l4'] = current_entity_set

		# entity data
		current_entity_info = {}
		self.layer_entity_info_sets['l4'] = current_entity_info

		# layer attributes
		self.layer_attributes['l4'] = {}
		################### end layer_info l4 ###################

