import datetime
import pprint

POINT = 1
EDGE = 2
FACE = 3
VOLUME = 4

def GetGeometryType(eltype):
    if eltype == 15:
        return POINT
    elif eltype in [1, 8, 26, 27, 28]:
        return EDGE
    elif eltype in [2, 3, 9, 10, 16, 20, 21, 22, 23, 24, 25]:
        return FACE
    elif eltype in [4, 5, 6, 7, 11, 12, 13, 14, 17, 18, 19, 29, 30, 31, 92, 93]:
        return VOLUME
    else:
        raise Exception("Unknown gmsh element type " + str(eltype))

def GetString(gtype):
    if gtype == POINT:
        return "POINT"
    elif gtype == EDGE:
        return "EDGE"
    elif gtype == FACE:
        return "FACE"
    elif gtype == VOLUME:
        return "VOLUME"
    else:
        return "Undefined"

def PostFix(dim, eltype):
    if dim == 2:
        pf = "2D"
    elif dim == 3:
        pf = "3D"
    if eltype in [1]:
        nn = 2
    elif eltype in [2, 8]:
        nn = 3
    elif eltype in [3, 4, 26]:
        nn = 4
    elif eltype in [6, 9, 28]:
        nn = 6
    elif eltype in [5, 16]:
        nn = 8
    elif eltype in [10, 20]:
        nn = 9
    elif eltype in [11, 21]:
        nn = 10
    elif eltype in [17, 29]:
        nn = 20
    elif eltype in [12]:
        nn = 27
    return pf + str(nn) + "N"

# node numbering conversion Gmsh-Kratos
def Permutation(eltype):
    if eltype == 1: # line2
        return [0, 1]
    elif eltype == 2: # tri3
        return [0, 1, 2]
    elif eltype == 3: # quad4
        return [0, 1, 2, 3]
    elif eltype == 4: # tet4
        return [0, 1, 2, 3]
    elif eltype == 5: # hex8
        return [0, 1, 2, 3, 4, 5, 6, 7]
    elif eltype == 8: # line3
        return [0, 1, 2]
    elif eltype == 9: # tri6
        return [0, 1, 2, 3, 4, 5]
    elif eltype == 10: # quad9
        return [0, 1, 2, 3, 4, 5, 6, 7, 8]
    elif eltype == 11: # tet10
        return [0, 1, 2, 3, 4, 5, 6, 7, 9, 8]
    elif eltype == 16: # quad8
        return [0, 1, 2, 3, 4, 5, 6, 7]
    elif eltype == 15: # point
        return [0]
    else:
        raise Exception("Unimplemented node numbering conversion for element type " + str(eltype))

class GmshMDPAWriter:
    """Class to write mdpa file from gmsh"""
    def __init__(self, dim, mesh):
        # save the gmsh mesh. mesh should be from gmshparse, e.g.
        # mesh = gmshparser.parse(gmsh_filename)
        self.mesh = mesh

        self.dim = dim
        self.element_tag_dict = {}
        self.element_gtype_dict = {}
        self.condition_tag_dict = {}
        self.condition_gtype_dict = {}

        self.layer_tag_dict = {}
        self.layer_gtype_dict = {}
        self.node_group_tag_dict = {}
        self.node_group_gtype_dict = {}
        self.element_tag_new_tag_dict = {}
        self.condition_tag_new_tag_dict = {}

    def Print(self):
        print("Summary of the existing tag in the mesh:")
        for entity in self.mesh.get_element_entities():
            eltype = entity.get_element_type()
            eltag = entity.get_tag()
            print("  Element tag: " + str(eltag) + ", type: " + str(eltype) + ", geometry type: " + str(GetString(GetGeometryType(eltype))))
        print("End of Summary")

    def SetElementTag(self, tag, gtype, element_name):
        self.element_gtype_dict[tag] = gtype
        self.element_tag_dict[tag] = element_name
        self.element_tag_new_tag_dict[tag] = tag

    def SetElementTagChange(self, tag, new_tag, gtype, element_name):
        self.element_gtype_dict[tag] = gtype
        self.element_tag_dict[tag] = element_name
        self.element_tag_new_tag_dict[tag] = new_tag

    def SetConditionTag(self, tag, gtype, condition_name):
        self.condition_gtype_dict[tag] = gtype
        self.condition_tag_dict[tag] = condition_name
        self.condition_tag_new_tag_dict[tag] = tag

    def SetConditionTagChange(self, tag, new_tag, gtype, condition_name):
        self.condition_gtype_dict[tag] = gtype
        self.condition_tag_dict[tag] = condition_name
        self.condition_tag_new_tag_dict[tag] = new_tag

    def SetLayerTag(self, tag, gtype, layer_name):
        self.layer_gtype_dict[tag] = gtype
        self.layer_tag_dict[tag] = layer_name

    def SetNodeGroupTag(self, tag, gtype, layer_name):
        self.node_group_gtype_dict[tag] = gtype
        self.node_group_tag_dict[tag] = layer_name

    def WriteMDPA(self, filename):
        ifile = open(filename + ".mdpa", "w")

        self.MDPA_Header(ifile)
        self.MDPA_Data(ifile)
        self.MDPA_Properties(ifile)
        self.MDPA_Nodes(ifile)
        self.MDPA_Elements(ifile)
        self.MDPA_Conditions(ifile)

        print("Writing MDPA to " + str(filename) + ".mdpa completed")
        ifile.close()

    def MDPA_Header(self, ifile):
        td = datetime.datetime.today()
        ifile.write("//KRATOS analysis data file\n")
        ifile.write("//(c) " + str(td.year) + " Hoang-Giang Bui, Ruhr-University Bochum\n")
        ifile.write("//This file is created on " + str(td.day) + "/" + str(td.month) + "/" + str(td.year % 100) + " " + str(td.hour) + ":" + str(td.minute) + ":" + str(td.second) + "\n")

    def MDPA_Data(self, ifile):
        ifile.write("\nBegin ModelPartData\n")
        ifile.write("End ModelPartData\n")

    def MDPA_Properties(self, ifile):
        tags = []
        for entity in self.mesh.get_element_entities():
            tags.append(entity.get_tag())
        for tag, new_tag in self.element_tag_new_tag_dict.items():
            tags.append(new_tag)
        for tag, new_tag in self.condition_tag_new_tag_dict.items():
            tags.append(new_tag)
        tags = list(set(tags))
        for tag in tags:
            ifile.write("\nBegin Properties " + str(tag) + "\n")
            ifile.write("End Properties\n")

    def MDPA_Nodes(self, ifile):
        ifile.write("\nBegin Nodes\n")
        for entity in self.mesh.get_node_entities():
            eltag = entity.get_tag()
            for node in entity.get_nodes():
                nid = node.get_tag()
                ncoords = node.get_coordinates()
                ifile.write(str(nid) + " " + str(ncoords[0]) + " " + str(ncoords[1]) + " " + str(ncoords[2]) + "\n")
            print("Wrote " + str(len(entity.get_nodes())) + " nodes from tag " + str(eltag))
        ifile.write("End Nodes\n")

    def MDPA_Elements(self, ifile):
        for entity in self.mesh.get_element_entities():
            eltype = entity.get_element_type()
            eltag = entity.get_tag()
            if eltag in self.element_tag_dict:
                gtype = self.element_gtype_dict[eltag]
                if GetGeometryType(eltype) == gtype:
                    elname = self.element_tag_dict[eltag]
                    conn_map = Permutation(eltype)
                    kel_name = elname + PostFix(self.dim, eltype)
                    ifile.write("\nBegin Elements " + kel_name + "\n")
                    for element in entity.get_elements():
                        elid = element.get_tag()
                        elcon = element.get_connectivity()
                        ifile.write(str(elid) + " " + str(self.element_tag_new_tag_dict[eltag]))
                        for i in range(0, len(elcon)):
                            ifile.write(" " + str(elcon[conn_map[i]]))
                        ifile.write("\n")
                    ifile.write("End Elements\n")
                    print("Wrote " + str(len(entity.get_elements())) + " elements of type " + kel_name + " of tag " + str(eltag) \
                        + ", Properties: " + str(self.element_tag_new_tag_dict[eltag]))

    def MDPA_Conditions(self, ifile):
        for entity in self.mesh.get_element_entities():
            eltype = entity.get_element_type()
            eltag = entity.get_tag()
            if eltag in self.condition_tag_dict:
                gtype = self.condition_gtype_dict[eltag]
                if GetGeometryType(eltype) == gtype:
                    coname = self.condition_tag_dict[eltag]
                    conn_map = Permutation(eltype)
                    kco_name = coname + PostFix(self.dim, eltype)
                    ifile.write("\nBegin Conditions " + kco_name + "\n")
                    for element in entity.get_elements():
                        elid = element.get_tag()
                        elcon = element.get_connectivity()
                        ifile.write(str(elid) + " " + str(self.condition_tag_new_tag_dict[eltag]))
                        for i in range(0, len(elcon)):
                            ifile.write(" " + str(elcon[conn_map[i]]))
                        ifile.write("\n")
                    ifile.write("End Conditions\n")
                    print("Wrote " + str(len(entity.get_elements())) + " conditions of type " + kco_name + " of tag " + str(eltag) \
                        + ", Properties: " + str(self.condition_tag_new_tag_dict[eltag]))

    def WriteLayers(self, filename):
        ifile = open(filename + "_layers.py", "w")

        layer_sets = {}
        layer_nodes_sets = {}

        for entity in self.mesh.get_element_entities():
            eltype = entity.get_element_type()
            eltag = entity.get_tag()
            if eltag in self.layer_tag_dict:
                gtype = self.layer_gtype_dict[eltag]
                if GetGeometryType(eltype) == gtype:
                    layer_name = self.layer_tag_dict[eltag]
                    if not (layer_name in layer_sets):
                        layer_sets[layer_name] = []
                    if not (layer_name in layer_nodes_sets):
                        layer_nodes_sets[layer_name] = []
                    for element in entity.get_elements():
                        elid = element.get_tag()
                        layer_sets[layer_name].append(elid)
                        elcon = element.get_connectivity()
                        layer_nodes_sets[layer_name].extend(elcon)

        for layer_name, layer_set in layer_sets.items():
            layer_sets[layer_name] = list(set(layer_set))

        ifile.write("def ReadLayerSets():\n")
        ifile.write("\tlayer_sets=" + str(layer_sets) + "\n")
        # ifile.write("\tlayer_sets=")
        # pprint.pprint(layer_sets, ifile)
        ifile.write("\treturn layer_sets\n")
        ifile.write("\n")

        # ifile.write("def ReadLayerSets():\n")
        # ifile.write("\tlayer_sets={}\n")
        # for entity in self.mesh.get_element_entities():
        #     eltype = entity.get_element_type()
        #     eltag = entity.get_tag()
        #     if eltag in self.layer_tag_dict:
        #         gtype = self.layer_gtype_dict[eltag]
        #         if GetGeometryType(eltype) == gtype:
        #             layer_name = self.layer_tag_dict[eltag]
        #             ifile.write("\tlayer_elements_list=[")
        #             cnt = 0
        #             for element in entity.get_elements():
        #                 elid = element.get_tag()
        #                 ifile.write("" + str(elid) + ",")
        #                 cnt = cnt+1
        #                 if cnt == 30:
        #                     ifile.write("\n\t")
        #                     cnt = 0
        #             ifile.write("]\n")
        #             ifile.write("\tlayer_sets['" + layer_name + "']=layer_elements_list\n")
        # ifile.write("\treturn layer_sets\n")
        # ifile.write("\n")

        for layer_name, layer_nodes_set in layer_nodes_sets.items():
            layer_nodes_sets[layer_name] = list(set(layer_nodes_set))

        ifile.write("def ReadLayerNodesSets():\n")
        ifile.write("\tlayer_nodes_sets=" + str(layer_nodes_sets) + "\n")
        # ifile.write("\tlayer_nodes_sets=")
        # pprint.pprint(tlayer_nodes_sets, ifile)
        ifile.write("\treturn layer_nodes_sets\n")
        ifile.write("\n")

        node_groups = {}

        for entity in self.mesh.get_element_entities():
            eltype = entity.get_element_type()
            eltag = entity.get_tag()
            if eltag in self.node_group_tag_dict:
                gtype = self.node_group_gtype_dict[eltag]
                if GetGeometryType(eltype) == gtype:
                    layer_name = self.node_group_tag_dict[eltag]
                    if not (layer_name in node_groups):
                        node_groups[layer_name] = []
                    for element in entity.get_elements():
                        elcon = element.get_connectivity()
                        node_groups[layer_name].extend(elcon)

        for layer_name, nodes in node_groups.items():
            node_groups[layer_name] = list(set(nodes))

        ifile.write("def ReadNodeGroups():\n")
        ifile.write("\tnode_groups=" + str(node_groups) + "\n")
        # ifile.write("\tnode_groups=")
        # pprint.pprint(tnode_groups, ifile)
        ifile.write("\treturn node_groups\n")
        ifile.write("\n")

        ifile.write("def ReadTopSurfaceNodes():\n")
        ifile.write("\treturn []\n")
        ifile.write("\n")

        ifile.write("def ReadBoundaryNodes():\n")
        ifile.write("\treturn []\n")
        ifile.write("\n")

        print("Writing layers to " + str(filename) + "_layers.py completed")
        ifile.close()

### Backup: Element type in Gmsh:

## elementType is e.g.:

# 1
# 2-node line.

# 2
# 3-node triangle.

# 3
# 4-node quadrangle.

# 4
# 4-node tetrahedron.

# 5
# 8-node hexahedron.

# 6
# 6-node prism.

# 7
# 5-node pyramid.

# 8
# 3-node second order line (2 nodes associated with the vertices and 1 with the edge).

# 9
# 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).

# 10
# 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).

# 11
# 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).

# 12
# 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).

# 13
# 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).

# 14
# 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).

# 15
# 1-node point.

# 16
# 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).

# 17
# 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).

# 18
# 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).

# 19
# 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).

# 20
# 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)

# 21
# 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)

# 22
# 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)

# 23
# 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)

# 24
# 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)

# 25
# 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)

# 26
# 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)

# 27
# 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)

# 28
# 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)

# 29
# 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)

# 30
# 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)

# 31
# 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)

# 92
# 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)

# 93
# 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)


# 9.2 Node ordering
# Historically, Gmsh first supported linear elements (lines, triangles, quadrangles, tetrahedra, prisms and hexahedra).
# Then, support for second and some third order elements has been added.
# Below we distinguish such "low order elements", which are hardcoded (i.e. they are explicitly defined in the code),
# and general "high-order elements", that have been coded in a more general fashion, theoretically valid for any order.

# 9.2.1 Low order elements
# For all mesh and post-processing file formats, the reference elements are defined as follows.

# Line:                 Line3:          Line4:

#       v
#       ^
#       |
#       |
# 0-----+-----1 --> u   0----2----1     0---2---3---1

# Triangle:               Triangle6:          Triangle9/10:          Triangle12/15:

# v
# ^                                                                   2
# |                                                                   | \
# 2                       2                    2                      9   8
# |`\                     |`\                  | \                    |     \
# |  `\                   |  `\                7   6                 10 (14)  7
# |    `\                 5    `4              |     \                |         \
# |      `\               |      `\            8  (9)  5             11 (12) (13) 6
# |        `\             |        `\          |         \            |             \
# 0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1

# Quadrangle:            Quadrangle8:            Quadrangle9:

#       v
#       ^
#       |
# 3-----------2          3-----6-----2           3-----6-----2
# |     |     |          |           |           |           |
# |     |     |          |           |           |           |
# |     +---- | --> u    7           5           7     8     5
# |           |          |           |           |           |
# |           |          |           |           |           |
# 0-----------1          0-----4-----1           0-----4-----1

# Tetrahedron:                          Tetrahedron10:
# (Note that for tet10, the node 8 and 9 is swapped in GiD/Kratos)
#                    v
#                  .
#                ,/
#               /
#            2                                     2
#          ,/|`\                                 ,/|`\
#        ,/  |  `\                             ,/  |  `\
#      ,/    '.   `\                         ,6    '.   `5
#    ,/       |     `\                     ,/       8     `\
#  ,/         |       `\                 ,/         |       `\
# 0-----------'.--------1 --> u         0--------4--'.--------1
#  `\.         |      ,/                 `\.         |      ,/
#     `\.      |    ,/                      `\.      |    ,9
#        `\.   '. ,/                           `7.   '. ,/
#           `\. |/                                `\. |/
#              `3                                    `3
#                 `\.
#                    ` w
# Hexahedron:             Hexahedron20:          Hexahedron27:
# (Note that for hex20 and hex27 the node numbering is different to GiD/Kratos)
#        v
# 3----------2            3----13----2           3----13----2
# |\     ^   |\           |\         |\          |\         |\
# | \    |   | \          | 15       | 14        |15    24  | 14
# |  \   |   |  \         9  \       11 \        9  \ 20    11 \
# |   7------+---6        |   7----19+---6       |   7----19+---6
# |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
# 0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
#  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
#   \ |     \  \ |         10 |        12|        10 |  21    12|
#    \|      w  \|           \|         \|          \|         \|
#     4----------5            4----16----5           4----16----5

# Prism:                      Prism15:               Prism18:

#            w
#            ^
#            |
#            3                       3                      3
#          ,/|`\                   ,/|`\                  ,/|`\
#        ,/  |  `\               12  |  13              12  |  13
#      ,/    |    `\           ,/    |    `\          ,/    |    `\
#     4------+------5         4------14-----5        4------14-----5
#     |      |      |         |      8      |        |      8      |
#     |    ,/|`\    |         |      |      |        |    ,/|`\    |
#     |  ,/  |  `\  |         |      |      |        |  15  |  16  |
#     |,/    |    `\|         |      |      |        |,/    |    `\|
#    ,|      |      |\        10     |      11       10-----17-----11
#  ,/ |      0      | `\      |      0      |        |      0      |
# u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
#     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
#     |,/         `\|         |,/         `\|        |,/         `\|
#     1-------------2         1------9------2        1------9------2

# Pyramid:                     Pyramid13:                   Pyramid14:

#                4                            4                            4
#              ,/|\                         ,/|\                         ,/|\
#            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
#          ,/   | | \                   ,/   | | \                   ,/   | | \
#        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
#      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
#    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
#  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
# 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
#  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
#    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
#      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \
#        `\.'      `\     `\`          `\.'             `\`         `\.'             `\`
#           1----------------2            1--------8-------2           1--------8-------2
#                     `\

# 9.2.2 High-order elements
# The node ordering of a higher order (possibly curved) element is compatible with the numbering of low order element (it is a generalization). We number nodes in the following order:

# - the element principal or corner vertices
# - the internal nodes for each edge
# - the internal nodes for each face
# - the volume internal nodes
# The numbering for internal nodes is recursive, i.e. the numbering follows that of the nodes of an embedded edge/face/volume of lower order.
# The higher order nodes are assumed to be equispaced. Edges and faces are numbered following the lowest order template that generates
# a single high-order on this edge/face. Furthermore, an edge is oriented from the node with the lowest to the highest index.
# The orientation of a face is such that the computed normal points outward; the starting point is the node with the lowest index.
