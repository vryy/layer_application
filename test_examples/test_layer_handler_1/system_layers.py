# geometry layers data file for isogeometric application
# (c) 2014 Hoang Giang Bui, Ruhr-University Bochum

class LayerProvider():
	def __init__(self):
		self.layer_list = []
		self.layer_attributes = {}
		self.layer_nodes_sets = {}
		self.layer_entities_sets = {}
		self.layer_entity_info_sets = {}

		## begin layer_info for layer sphere
		self.layer_list.append('sphere')

		# nodal info
		current_nodal_set = {}
		current_nodal_set[1] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[2] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[3] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[4] = [2.000000, 0.000000, 0.000000]
		current_nodal_set[5] = [2.000000, 0.000000, 0.000000]
		current_nodal_set[6] = [2.000000, 0.000000, 0.000000]
		current_nodal_set[7] = [1.000000, 1.000000, 0.000000]
		current_nodal_set[8] = [1.000000, 1.000000, 1.000000]
		current_nodal_set[9] = [1.000000, 0.000000, 1.000000]
		current_nodal_set[10] = [2.000000, 2.000000, 0.000000]
		current_nodal_set[11] = [2.000000, 2.000000, 2.000000]
		current_nodal_set[12] = [2.000000, 0.000000, 2.000000]
		current_nodal_set[13] = [0.000000, 1.000000, -0.000000]
		current_nodal_set[14] = [-0.000000, 1.000000, 1.000000]
		current_nodal_set[15] = [-0.000000, 0.000000, 1.000000]
		current_nodal_set[16] = [0.000000, 2.000000, -0.000000]
		current_nodal_set[17] = [-0.000000, 2.000000, 2.000000]
		current_nodal_set[18] = [-0.000000, 0.000000, 2.000000]
		self.layer_nodes_sets['sphere'] = current_nodal_set

		# entity connectivities
		current_entity_set = {}
		current_entity_set[1] = [1, 7, 13, 4, 10, 16, 2, 8, 14, 5, 11, 17, 3, 9, 15, 6, 12, 18]
		self.layer_entities_sets['sphere'] = current_entity_set

		# entity data
		current_entity_info = {}
		temp = {}
		temp[1] = [1.000000, 0.707107, 1.000000, 1.000000, 0.707107, 1.000000, 0.707107, 0.500000, 0.707107, 0.707107, 0.500000, 0.707107, 1.000000, 0.707107, 1.000000, 1.000000, 0.707107, 1.000000]
		current_entity_info['NURBS_WEIGHT'] = temp

		temp = {}
		temp[1] = [[ 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19],[ 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 0.000000]]
		current_entity_info['EXTRACTION_OPERATOR'] = temp

		temp = {}
		temp[1] = 2
		current_entity_info['NURBS_DEGREE_1'] = temp

		temp = {}
		temp[1] = 1
		current_entity_info['NURBS_DEGREE_2'] = temp

		temp = {}
		temp[1] = 2
		current_entity_info['NURBS_DEGREE_3'] = temp

		temp = {}
		temp[1] = 1
		current_entity_info['NUM_DIVISION_1'] = temp

		temp = {}
		temp[1] = 1
		current_entity_info['NUM_DIVISION_2'] = temp

		temp = {}
		temp[1] = 1
		current_entity_info['NUM_DIVISION_3'] = temp

		self.layer_entity_info_sets['sphere'] = current_entity_info

		# layer attributes
		self.layer_attributes['sphere'] = {}
		################### end layer_info sphere ###################

		## begin layer_info for layer surface
		self.layer_list.append('surface')

		# nodal info
		current_nodal_set = {}
		current_nodal_set[1] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[2] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[3] = [1.000000, 0.000000, 0.000000]
		current_nodal_set[4] = [1.000000, 1.000000, 0.000000]
		current_nodal_set[5] = [1.000000, 1.000000, 1.000000]
		current_nodal_set[6] = [1.000000, 0.000000, 1.000000]
		current_nodal_set[7] = [0.000000, 1.000000, -0.000000]
		current_nodal_set[8] = [-0.000000, 1.000000, 1.000000]
		current_nodal_set[9] = [-0.000000, 0.000000, 1.000000]
		self.layer_nodes_sets['surface'] = current_nodal_set

		# entity connectivities
		current_entity_set = {}
		current_entity_set[1] = [1, 4, 7, 2, 5, 8, 3, 6, 9]
		self.layer_entities_sets['surface'] = current_entity_set

		# entity data
		current_entity_info = {}
		temp = {}
		temp[1] = [1.000000, 0.707107, 1.000000, 0.707107, 0.500000, 0.707107, 1.000000, 0.707107, 1.000000]
		current_entity_info['NURBS_WEIGHT'] = temp

		temp = {}
		temp[1] = [[ 10, 10, 10, 10, 10, 10, 10, 10, 10, 10],[ 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 0.000000]]
		current_entity_info['EXTRACTION_OPERATOR'] = temp

		temp = {}
		temp[1] = 2
		current_entity_info['NURBS_DEGREE_1'] = temp

		temp = {}
		temp[1] = 2
		current_entity_info['NURBS_DEGREE_2'] = temp

		temp = {}
		temp[1] = 1
		current_entity_info['NUM_DIVISION_1'] = temp

		temp = {}
		temp[1] = 1
		current_entity_info['NUM_DIVISION_2'] = temp

		self.layer_entity_info_sets['surface'] = current_entity_info

		# layer attributes
		self.layer_attributes['surface'] = {}
		################### end layer_info surface ###################

	def CollapseLayer(self, layer_list):
		import math
		tol = 1.0e-6

		# collapse each layer specified in layer list
		for str_layer in  layer_list:
			node_map = {}
			non_repetitive_nodes = []
			# iterate all the nodes to check for node repetition
			for i_node in self.layer_nodes_sets[str_layer]:
				n = self.layer_nodes_sets[str_layer][i_node]
				for i_other_node in self.layer_nodes_sets[str_layer]:
					n1 = self.layer_nodes_sets[str_layer][i_other_node]
					d = math.sqrt(math.pow(n[0] - n1[0], 2) + math.pow(n[1] - n1[1], 2) + math.pow(n[2] - n1[2], 2))
					if d < tol:
						node_map[i_node] = i_other_node
						if i_other_node == i_node:
							non_repetitive_nodes.append(i_node)
						break

			# reform the layer nodal set
			new_nodes_set = {}
			node_map_nonrepetitive = {}
			cnt = 1
			for i_node in non_repetitive_nodes:
				new_nodes_set[cnt] = self.layer_nodes_sets[str_layer][i_node]
				node_map_nonrepetitive[i_node] = cnt
				cnt = cnt + 1
			self.layer_nodes_sets[str_layer] = new_nodes_set

			# reform the entity connectivities
			for i_entity in self.layer_entities_sets[str_layer]:
				new_entity = []
				for i_node in self.layer_entities_sets[str_layer][i_entity]:
					new_entity.append(node_map_nonrepetitive[node_map[i_node]])
				self.layer_entities_sets[str_layer][i_entity] = new_entity

	def GlueLayer(self, str_source_layer, str_target_layer):
		import math
		tol = 1.0e-6

		# signal that the str_source_layer has been glued to str_target_layer
		self.layer_attributes[str_source_layer]['master'] = str_target_layer
		
		# first check the largest nodal id in target_layer
		i_max_node = max(self.layer_nodes_sets[str_target_layer])

		# iterate all the nodes in source_layer to find the repetitive nodes
		node_map = {}
		for i_node in self.layer_nodes_sets[str_source_layer]:
			n = self.layer_nodes_sets[str_source_layer][i_node]
			for i_other_node in self.layer_nodes_sets[str_target_layer]:
				n1 = self.layer_nodes_sets[str_target_layer][i_other_node]
				d = math.sqrt(math.pow(n[0] - n1[0], 2) + math.pow(n[1] - n1[1], 2) + math.pow(n[2] - n1[2], 2))
				if d < tol:
					node_map[i_node] = i_other_node
					break
		print "node_map:", node_map

		# reform the nodal set of source_layer
		new_nodes_set = {}
		cnt = i_max_node + 1
		for i_node in self.layer_nodes_sets[str_source_layer]:
			if i_node not in node_map:
				new_nodes_set[cnt] = self.layer_nodes_sets[str_source_layer][i_node]
				node_map[i_node] = cnt
				cnt = cnt + 1
		self.layer_nodes_sets[str_source_layer] = new_nodes_set

		# reform the entity set of source_layer
		for i_entity in self.layer_entities_sets[str_source_layer]:
			new_entity = []
			for i_node in self.layer_entities_sets[str_source_layer][i_entity]:
				new_entity.append(node_map[i_node])
			self.layer_entities_sets[str_source_layer][i_entity] = new_entity

	def RenumberAll(self):
		cnt = 1
		for str_layer in self.layer_list:
			if 'master' not in self.layer_attributes[str_layer]:
				# renumber all nodes in current layer
				new_nodes_set = {}
				node_map = {}
				for i_node in self.layer_nodes_sets[str_layer]:
					new_nodes_set[cnt] = self.layer_nodes_sets[str_layer][i_node]
					node_map[i_node] = cnt
					cnt = cnt + 1
				self.layer_nodes_sets[str_layer] = new_nodes_set
				print "node_map:", node_map

				# reform the entity connectivities in current layer
				for i_entity in self.layer_entities_sets[str_layer]:
					new_entity = []
					for i_node in self.layer_entities_sets[str_layer][i_entity]:
						new_entity.append(node_map[i_node])
					self.layer_entities_sets[str_layer][i_entity] = new_entity

#import pprint
#layers = LayerProvider()
#pp = pprint.PrettyPrinter(indent=4,depth=6)
#pp.pprint(layers.layer_nodes_sets)
#pp.pprint(layers.layer_entities_sets)
#pp.pprint(layers.layer_entity_info_sets)
