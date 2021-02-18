class FEBioWriter:
    """Class to write feb file from Kratos model_part"""
    def __init__(self):
        self._Geometry_Tag_ = "Geometry"
        self._Nodes_Tag_ = "Nodes"
        self._Elements_Tag_ = "Elements"
        self._NodeSet_Tag_ = "NodeSet"
        self._Surface_Tag_ = "Surface"
        self._Material_Tag_ = "Material"

        max_tags = 10
        self.tags = []
        for i in range(0, max_tags):
            self.tags.append("")

        self.gid2febio_hex20 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15]
        self.gid2febio_hex27 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15, 21, 22, 23, 24, 20, 25, 26]

    def BeginFebruary(self, filename, spec_number, module_type):
        self.level = 0
        self.ifile = open(filename + ".feb", "w")
        self.ifile.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        self.ifile.write('<febio_spec version="' + str(spec_number) + '">\n')
        self.ifile.write('\t<Module type="' + module_type + '"/>\n')

    def BeginGeometry(self):
        self.tags[self.level] = self._Geometry_Tag_
        current_tag = self.tags[self.level]
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        self.ifile.write(tab + '<' + current_tag + '>\n')

    def BeginNodes(self, set_name):
        self.tags[self.level] = self._Nodes_Tag_
        current_tag = self.tags[self.level]
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        self.ifile.write(tab + '<' + current_tag + ' name="' + set_name + '">\n')

    def WriteNodes(self, nodes):
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        for node in nodes:
            self.ifile.write(tab + '<node id="' + str(node.Id) + '">')
            self.ifile.write(str(node.X0) + ', ' + str(node.Y0) + ', ' + str(node.Z0))
            self.ifile.write('</node>\n')
        self.level = self.level - 1

    def WriteNodesFromList(self, list_nodes):
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        for node_id, coord in list_nodes.iteritems():
            self.ifile.write(tab + '<node id="' + str(node_id) + '">')
            self.ifile.write(str(coord[0]) + ', ' + str(coord[1]) + ', ' + str(coord[2]))
            self.ifile.write('</node>\n')
        self.level = self.level - 1

    def EndNodes(self):
        tab = self.GetTabLevel(self.level)
        self.level = self.level - 1
        current_tag = self.tags[self.level]
        if current_tag != self._Nodes_Tag_:
            raise Exception("ERROR: Begin" + current_tag + " is not yet called\n")
        self.ifile.write(tab + '</' + current_tag + '>\n')

    def BeginElements(self, type, mat_id, set_name):
        self.tags[self.level] = self._Elements_Tag_
        current_tag = self.tags[self.level]
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        self.current_element_type = type
        self.ifile.write(tab + '<' + current_tag + ' type="' + type + '" mat="' + str(mat_id) + '"' + ' name="' + set_name + '">\n')

    def WriteElements(self, elements):
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        for elem in elements:
            list_nodes = []
            for node in elem.GetNodes():
                list_nodes.append(node.Id)
            if self.current_element_type == "hex27":
                new_list_nodes = []
                for i in range(0, 27):
                    new_list_nodes.append(list_nodes[self.gid2febio_hex27[i]])
                list_nodes = new_list_nodes
            if self.current_element_type == "hex20":
                new_list_nodes = []
                for i in range(0, 20):
                    new_list_nodes.append(list_nodes[self.gid2febio_hex20[i]])
                list_nodes = new_list_nodes
            self.ifile.write(tab + '<elem id="' + str(elem.Id) + '">')
            first = True
            for node_id in list_nodes:
                if first:
                    self.ifile.write(str(node_id))
                    first = False
                else:
                    self.ifile.write(', ' + str(node_id))
            self.ifile.write('</elem>\n')
        self.level = self.level - 1

    def EndElements(self):
        tab = self.GetTabLevel(self.level)
        self.level = self.level - 1
        current_tag = self.tags[self.level]
        if current_tag != self._Elements_Tag_:
            raise Exception("ERROR: Begin" + current_tag + " is not yet called\n")
        self.ifile.write(tab + '</' + current_tag + '>\n')

    def BeginNodeSet(self, set_name):
        self.tags[self.level] = self._NodeSet_Tag_
        current_tag = self.tags[self.level]
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        self.ifile.write(tab + '<' + current_tag + ' name="' + set_name + '">\n')

    def WriteNodeSet(self, list_nodes):
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        for node_id in list_nodes:
            self.ifile.write(tab + '<node id="' + str(node_id) + '"/>\n')
        self.level = self.level - 1

    def EndNodeSet(self):
        tab = self.GetTabLevel(self.level)
        self.level = self.level - 1
        current_tag = self.tags[self.level]
        if current_tag != self._NodeSet_Tag_:
            raise Exception("ERROR: Begin" + current_tag + " is not yet called\n")
        self.ifile.write(tab + '</' + current_tag + '>\n')

    def BeginSurface(self, surface_name, surface_type):
        self.tags[self.level] = self._Surface_Tag_
        current_tag = self.tags[self.level]
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        self.ifile.write(tab + '<' + current_tag + ' name="' + surface_name + '">\n')
        self.surface_type = surface_type

    def WriteSurface(self, conds):
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        for cond in conds:
            self.ifile.write(tab + '<' + self.surface_type + ' id="' + str(cond.Id) + '">')
            first = True
            for node in cond.GetNodes():
                if first:
                    self.ifile.write(str(node.Id))
                    first = False
                else:
                    self.ifile.write(', ' + str(node.Id))
            self.ifile.write('</' + self.surface_type + '>\n')
        self.level = self.level - 1

    def EndSurface(self):
        tab = self.GetTabLevel(self.level)
        self.level = self.level - 1
        current_tag = self.tags[self.level]
        if current_tag != self._Surface_Tag_:
            raise Exception("ERROR: Begin" + current_tag + " is not yet called\n")
        self.ifile.write(tab + '</' + current_tag + '>\n')

    def EndGeometry(self):
        tab = self.GetTabLevel(self.level)
        self.level = self.level - 1
        current_tag = self.tags[self.level]
        if current_tag != self._Geometry_Tag_:
            raise Exception("ERROR: Begin" + current_tag + " is not yet called\n")
        self.ifile.write(tab + '</' + current_tag + '>\n')

    def BeginMaterial(self):
        self.tags[self.level] = self._Material_Tag_
        current_tag = self.tags[self.level]
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        self.ifile.write(tab + '<' + current_tag + '>\n')

    # write a single material
    def WriteMaterial(self, mat_info):
        self.level = self.level + 1
        tab = self.GetTabLevel(self.level)
        tab1 = self.GetTabLevel(self.level+1)
        self.ifile.write(tab + '<material id="' + str(mat_info['id']) + '" name="' + str(mat_info['name']) + '" type="' + str(mat_info['type']) + '">\n')
        for tag, value in mat_info['data'].iteritems():
            self.ifile.write(tab1 + '<' + tag + '>' + str(value) + '</' + tag + '>\n')
        self.ifile.write(tab + '</material>\n')
        self.level = self.level - 1

    def EndMaterial(self):
        tab = self.GetTabLevel(self.level)
        self.level = self.level - 1
        current_tag = self.tags[self.level]
        if current_tag != self._Material_Tag_:
            raise Exception("ERROR: Begin" + current_tag + " is not yet called\n")
        self.ifile.write(tab + '</' + current_tag + '>\n')

    def EndFebruary(self):
        self.ifile.write('</febio_spec>\n')
        self.ifile.close()

    def GetTabLevel(self, nlevel):
        tab = ""
        for i in range(0, nlevel):
            tab = tab + '\t'
        return tab
