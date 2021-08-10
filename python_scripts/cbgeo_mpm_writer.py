from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *

class CbGeoMPMWriter:
    """Class to write feb file from Kratos model_part"""
    def __init__(self, dim):
        self.dim = dim
        self.offset = 1
        self.shift = 0

    def WriteMesh(self, model_part, filename):
        ifile = open(filename, "w")
        ifile.write(str(len(model_part.Nodes)) + " " + str(len(model_part.Elements)) + "\n")
        for node in model_part.Nodes:
            if self.dim == 2:
                ifile.write(str(node.X0) + " " + str(node.Y0) + "\n")
            elif self.dim == 3:
                ifile.write(str(node.X0) + " " + str(node.Y0) + " " + str(node.Z0) + "\n")
        for elem in model_part.Elements:
            nodes = []
            for node in elem.GetNodes():
                nodes.append(node.Id)
            for i in range(0, len(nodes)):
                if i > 0:
                    ifile.write(" ")
                ifile.write(str(nodes[(i+self.shift)%len(nodes)]-self.offset))
            ifile.write("\n")
        ifile.close()
        print("Mesh is written to " + filename + " successfully")

    def WriteParticles(self, model_part, filename):
        ifile = open(filename, "w")
        particles = []
        for elem in model_part.Elements:
            if elem.Is(ACTIVE):
                ipoints = elem.CalculateOnIntegrationPoints(INTEGRATION_POINT_GLOBAL, model_part.ProcessInfo)
                # print(ipoints)
                for p in ipoints:
                    particles.append(p)
        ifile.write(str(len(particles)) + "\n")
        for p in particles:
            ifile.write(str(p[0]))
            for i in range(1, self.dim):
                ifile.write(" " + str(p[i]))
            ifile.write("\n")
        ifile.close()
        print("Particles are written to " + filename + " successfully")

    def WriteParticlesVolumes(self, model_part, filename):
        ifile = open(filename, "w")
        for elem in model_part.Elements:
            if elem.Is(ACTIVE):
                jacobians = elem.CalculateOnIntegrationPoints(JACOBIAN_0, model_part.ProcessInfo)
                weights = elem.CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, model_part.ProcessInfo)
                for i in range(0, len(jacobians)):
                    ifile.write(str(jacobians[i][0]*weights[i][0]) + "\n")
        ifile.close()
        print("Particles volumes are written to " + filename + " successfully")

    def WriteEntitySets(self, layer_particle_sets, layer_node_sets, filename):
        ifile = open(filename, "w")
        ifile.write('{\n')
        ifile.write('  "particle_sets": [\n')
        particle_ini = True
        for index in layer_particle_sets.keys():
            if not particle_ini:
                ifile.write('    ,\n')
            else:
                particle_ini = False
            ifile.write('    {\n')
            ifile.write('      "id" : ' + str(index) + ',\n')
            ifile.write('      "set" : [\n')
            ifile.write('        ')
            index_ini = True
            for n in layer_particle_sets[index]:
                if not index_ini:
                    ifile.write(', ' + str(n-self.offset))
                else:
                    ifile.write(str(n-self.offset))
                    index_ini = False
            ifile.write('\n')
            ifile.write('      ]\n')
            ifile.write('    }\n')
        ifile.write('  ],\n')
        ifile.write('  "node_sets": [\n')
        node_ini = True
        for index in layer_node_sets.keys():
            if not node_ini:
                ifile.write('    ,\n')
            else:
                node_ini = False
            ifile.write('    {\n')
            ifile.write('      "id" : ' + str(index) + ',\n')
            ifile.write('      "set" : [\n')
            ifile.write('        ')
            index_ini = True
            for n in layer_node_sets[index]:
                if not index_ini:
                    ifile.write(', ' + str(n-self.offset))
                else:
                    ifile.write(str(n-self.offset))
                    index_ini = False
            ifile.write('\n')
            ifile.write('      ]\n')
            ifile.write('    }\n')
        ifile.write('  ]\n')
        ifile.write("}\n")
        ifile.close()
        print("Entity sets are written to " + filename + " successfully")

