#include <cstdlib>

///Flags for mesh writing
enum VTK_PostMode {VTK_PostAscii, VTK_PostBinary};

FILE* VTK_fOpenPostResultFile(const char* filename, VTK_PostMode mode)
{
    FILE* file = fopen(filename, "w");

    if(mode == VTK_PostAscii)
    {
        fprintf(file, "<?xml version=\"1.0\"?>\n");
        fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
        fprintf(file, "  <UnstructuredGrid>\n");
    }
    else if(mode == VTK_PostBinary)
    {
        fprintf(file, "<?xml version=\"1.0\"?>\n");
        fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n");
        fprintf(file, "  <UnstructuredGrid>\n");
    }

    return file;
}

void VTK_fClosePostResultFile( FILE* file )
{
    fprintf(file, "  </UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");
    fclose(file);
}

void VTK_fBeginMesh(FILE* file, const char* name, unsigned int num_points, unsigned int num_cells)
{
    fprintf(file, "    <Piece Name=\"%s\" NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", name, num_points, num_cells);
}

void VTK_fEndMesh(FILE* file)
{
    fprintf(file, "    </Piece>\n");
}

void VTK_fBeginCoordinates(FILE* file)
{
    fprintf(file, "      <Points>\n");
    fprintf(file, "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
}

void VTK_fEndCoordinates(FILE* file)
{
    fprintf(file, "        </DataArray>\n");
    fprintf(file, "      </Points>\n");
}

void VTK_fWriteCoordinates(FILE* file, unsigned int id, double x, double y, double z)
{
    fprintf(file, "          %f %f %f\n", x, y, z);
}

void VTK_fBeginElements(FILE* file)
{
    fprintf(file, "      <Cells>\n");
}

void VTK_fEndElements(FILE* file)
{
    fprintf(file, "      </Cells>\n");
}

void VTK_fBeginElementsConnectivity(FILE* file)
{
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"connectivity\">\n");
}

void VTK_fEndElementsConnectivity(FILE* file)
{
    fprintf(file, "        </DataArray>\n");
}

void VTK_fBeginElementsOffsets(FILE* file)
{
    fprintf(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
}

void VTK_fEndElementsOffsets(FILE* file)
{
    fprintf(file, "</DataArray>\n");
}

void VTK_fBeginElementsTypes(FILE* file)
{
    fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
}

void VTK_fEndElementsTypes(FILE* file)
{
    fprintf(file, "</DataArray>\n");
}

void VTK_fWriteElementConnectivity(FILE* file, unsigned int id, int* nodes_id, unsigned int nodes_size)
{
    unsigned int i;
    fprintf(file, "         ");
    for(i = 0; i < nodes_size; ++i)
        fprintf(file, " %d", nodes_id[i]);
    fprintf(file, "\n");
}

void VTK_fWriteElementOffset(FILE* file, unsigned int offset)
{
    fprintf(file, " %d", offset);
}

void VTK_fWriteElementType(FILE* file, unsigned int type)
{
    fprintf(file, " %d", type);
}

void VTK_fBeginResult(FILE* file, const char* result_name, unsigned int component_size, VTK_PostMode mode)
{
    if(mode == VTK_PostAscii)
    {
        fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", result_name, component_size);
    }
    else if(mode == VTK_PostBinary)
    {
        fprintf(file, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"binary\">\n", result_name, component_size);
    }
}

void VTK_fEndResult(FILE* file)
{
    fprintf(file, "        </DataArray>\n");
}

void VTK_fWriteScalar(FILE* file, unsigned int id, double value)
{
    fprintf(file, "          %f\n", value);
}

void VTK_fWriteVector3(FILE* file, unsigned int id, double value1, double value2, double value3)
{
    fprintf(file, "          %f %f %f\n", value1, value2, value3);
}

void VTK_fWriteVector(FILE* file, unsigned int id, unsigned int size, double* values)
{
    fprintf(file, "         ");
    for(unsigned int i = 0; i < size; ++i)
        fprintf(file, " %f", values[i]);
    fprintf(file, "\n");
}


