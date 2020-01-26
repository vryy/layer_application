#if !defined(KRATOS_VTK_H_INCLUDED)
#define  KRATOS_VTK_H_INCLUDED

#include <cstdlib>

#if defined(_USE_ZLIB) && defined(_USE_LIBB64)
#include "external_libraries/libb64/libb64.h"
#include "external_libraries/zlib/zlib.h"
//#include "zconf.h"
#endif

///Flags for mesh writing
enum VTK_PostMode {VTK_PostAscii, VTK_PostBinary};
enum VTK_PostFileFormat {VTK_PostVTU, VTK_PostVTM};

// reference: vtkCellType.h
enum VTK_ElementType {
    VTK_Empty_Cell = 0,
    VTK_Vertex = 1,
    VTK_Poly_Vertex = 2,
    VTK_Line = 3,
    VTK_Poly_Line = 4,
    VTK_Triangle = 5,
    VTK_Triangle_Strip = 6,
    VTK_Polygon = 7,
    VTK_Pixel = 8,
    VTK_Quad = 9,
    VTK_Tetra = 10,
    VTK_Voxel = 11,
    VTK_Hexahedron = 12,
    VTK_Wedge = 13,
    VTK_Pyramid = 14,
    VTK_Pentagonal_Prism = 15,
    VTK_Hexagonal_Prism = 16,
    VTK_Quadratic_Edge = 21,
    VTK_Quadratic_Triangle = 22,
    VTK_Quadratic_Quad = 23,
    VTK_Quadratic_Polygon = 36,
    VTK_Quadratic_Tetra = 24,
    VTK_Quadratic_Hexahedron = 25,
    VTK_Quadratic_Wedge = 26,
    VTK_Quadratic_Pyramid = 27,
    VTK_Biquadratic_Quad = 28,
    VTK_Triquadratic_Hexahedron = 29,
    VTK_Quadratic_Linear_Quad = 30,
    VTK_Quadratic_Linear_Wedge = 31,
    VTK_Biquadratic_Quadratic_Wedge = 32,
    VTK_Biquadratic_Quadratic_Hexahedron = 33,
    VTK_Biquadratic_Triangle = 34,
    VTK_Cubix_Line = 35,
    VTK_Convex_Point_Set = 41,
    VTK_Polyhedron = 42,
    VTK_Parametric_Curve = 51,
    VTK_Parametric_Surface = 52,
    VTK_Parametric_Tri_Surface = 53,
    VTK_Parametric_Quad_Surface = 54,
    VTK_Parametric_Tetra_Region = 55,
    VTK_Parametric_Hex_Region = 56,
    VTK_Higher_Order_Edge = 60,
    VTK_Higher_Order_Triangle = 61,
    VTK_Higher_Order_Quad = 62,
    VTK_Higher_Order_Polygon = 63,
    VTK_Higher_Order_Tetrahedron = 64,
    VTK_Higher_Order_Wedge = 65,
    VTK_Higher_Order_Pyramid = 66,
    VTK_Higher_Order_Hexahedron = 67,
    VTK_Number_Of_Cell_Types
};

FILE* VTK_fOpenPostVTUFile(const char* filename, VTK_PostMode mode)
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

void VTK_fClosePostVTUFile( FILE* file )
{
    fprintf(file, "  </UnstructuredGrid>\n");
    fprintf(file, "</VTKFile>\n");
    fclose(file);
}

FILE* VTK_fOpenPostVTMFile(const char* filename, VTK_PostMode mode)
{
    FILE* file = fopen(filename, "w");

    if(mode == VTK_PostAscii)
    {
        fprintf(file, "<?xml version=\"1.0\"?>\n");
        fprintf(file, "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\">\n");
        fprintf(file, "  <vtkMultiBlockDataSet>\n");
    }
    else if(mode == VTK_PostBinary)
    {
        fprintf(file, "<?xml version=\"1.0\"?>\n");
        fprintf(file, "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n");
        fprintf(file, "  <vtkMultiBlockDataSet>\n");
    }

    return file;
}

void VTK_fClosePostVTMFile( FILE* file )
{
    fprintf(file, "  </vtkMultiBlockDataSet>\n");
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

void VTK_fBeginCoordinates(FILE* file, VTK_PostMode mode)
{
    fprintf(file, "      <Points>\n");
    if (mode == VTK_PostAscii)
        fprintf(file, "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    else if (mode == VTK_PostBinary)
        fprintf(file, "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"binary\">\n");
}

void VTK_fEndCoordinates(FILE* file)
{
    fprintf(file, "        </DataArray>\n");
    fprintf(file, "      </Points>\n");
}

void VTK_fWriteCoordinates(FILE* file, unsigned int id, double x, double y, double z)
{
    fprintf(file, "          %e %e %e\n", x, y, z);
}

void VTK_fBeginElements(FILE* file)
{
    fprintf(file, "      <Cells>\n");
}

void VTK_fEndElements(FILE* file)
{
    fprintf(file, "      </Cells>\n");
}

void VTK_fBeginElementsConnectivity(FILE* file, VTK_PostMode mode)
{
    if (mode == VTK_PostAscii)
        fprintf(file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    else if (mode == VTK_PostBinary)
        fprintf(file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
}

void VTK_fEndElementsConnectivity(FILE* file)
{
    fprintf(file, "        </DataArray>\n");
}

void VTK_fBeginElementsOffsets(FILE* file, VTK_PostMode mode)
{
    if (mode == VTK_PostAscii)
        fprintf(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    else if (mode == VTK_PostBinary)
        fprintf(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
}

void VTK_fEndElementsOffsets(FILE* file)
{
    fprintf(file, "\n        </DataArray>\n");
}

void VTK_fBeginElementsTypes(FILE* file, VTK_PostMode mode)
{
    if (mode == VTK_PostAscii)
        fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    else if (mode == VTK_PostBinary)
        fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n");
}

void VTK_fEndElementsTypes(FILE* file)
{
    fprintf(file, "\n        </DataArray>\n");
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
    fprintf(file, "          %e\n", value);
}

void VTK_fWriteVector3(FILE* file, unsigned int id, double value1, double value2, double value3)
{
    fprintf(file, "          %e %e %e\n", value1, value2, value3);
}

void VTK_fWriteVector(FILE* file, unsigned int id, unsigned int size, double* values)
{
    fprintf(file, "         ");
    for(unsigned int i = 0; i < size; ++i)
        fprintf(file, " %e", values[i]);
    fprintf(file, "\n");
}

void VTK_fBeginDataset(FILE* file, const char* name, const char* filename, unsigned int index)
{
    fprintf(file, "    <DataSet index=\"%d\" name=\"%s\" file=\"%s\">\n", index, name, filename);
}

void VTK_fEndDataset(FILE* file)
{
    fprintf(file, "    </DataSet>\n");
}

#if defined(_USE_ZLIB) && defined(_USE_LIBB64)

void zlib_check(int val)
{
    if(val != Z_OK)
    {
        std::stringstream ss;
        ss << "Error with zlib, error code = " << val;
        throw std::runtime_error(ss.str());
    }
}

int vtk_write_compressed (FILE * vtkfile, char *numeric_data, size_t byte_length)
{
    int                 retval, fseek1, fseek2;
    size_t              iz;
    size_t              blocksize, lastsize;
    size_t              theblock, numregularblocks, numfullblocks;
    size_t              header_entries, header_size;
    size_t              code_length, base_length;
    long                header_pos, final_pos;
    char               *comp_data, *base_data;
    uint32_t           *compression_header;
    uLongf              comp_length;
    base64_encodestate  encode_state;

    /* compute block sizes */
    blocksize = (size_t) (1 << 15);       /* 32768 */
    lastsize = byte_length % blocksize;
    numregularblocks = byte_length / blocksize;
    numfullblocks = numregularblocks + (lastsize > 0 ? 1 : 0);
    header_entries = 3 + numfullblocks;
    header_size = header_entries * sizeof (uint32_t);

    /* allocate compression and base64 arrays */
    code_length = 2 * std::max (blocksize, header_size) + 4 + 1;
    comp_data = (char*) malloc(sizeof(char) * code_length);
    base_data = (char*) malloc(sizeof(char) * code_length);

    /* figure out the size of the header and write a dummy */
    compression_header = (uint32_t*) malloc(sizeof(uint32_t) * header_entries);
    compression_header[0] = (uint32_t) numfullblocks;
    compression_header[1] = (uint32_t) blocksize;
    compression_header[2] = (uint32_t)
    (lastsize > 0 || byte_length == 0 ? lastsize : blocksize);
    for (iz = 3; iz < header_entries; ++iz) {
        compression_header[iz] = 0;
    }
    base64_init_encodestate (&encode_state);
    base_length = base64_encode_block ((char *) compression_header,
                                     header_size, base_data, &encode_state);
    base_length +=
    base64_encode_blockend (base_data + base_length, &encode_state);
    assert (base_length < code_length);
    base_data[base_length] = '\0';
    header_pos = ftell (vtkfile);
    (void) fwrite (base_data, 1, base_length, vtkfile);

    /* write the regular data blocks */
    base64_init_encodestate (&encode_state);
    for (theblock = 0; theblock < numregularblocks; ++theblock) {
        comp_length = code_length;
        retval = compress2 ((Bytef *) comp_data, &comp_length,
                            (const Bytef *) (numeric_data + theblock * blocksize),
                            (uLong) blocksize, Z_BEST_COMPRESSION);
        zlib_check(retval);
        compression_header[3 + theblock] = comp_length;
        base_length = base64_encode_block (comp_data, comp_length,
                                           base_data, &encode_state);
        assert (base_length < code_length);
        base_data[base_length] = '\0';
        (void) fwrite (base_data, 1, base_length, vtkfile);
    }

    /* write odd-sized last block if necessary */
    if (lastsize > 0) {
        comp_length = code_length;
        retval = compress2 ((Bytef *) comp_data, &comp_length,
                            (const Bytef *) (numeric_data + theblock * blocksize),
                            (uLong) lastsize, Z_BEST_COMPRESSION);
        zlib_check(retval);
        compression_header[3 + theblock] = comp_length;
        base_length = base64_encode_block (comp_data, comp_length, base_data, &encode_state);
        assert (base_length < code_length);
        base_data[base_length] = '\0';
        (void) fwrite (base_data, 1, base_length, vtkfile);
    }

    /* write base64 end block */
    base_length = base64_encode_blockend (base_data, &encode_state);
    assert (base_length < code_length);
    base_data[base_length] = '\0';
    (void) fwrite (base_data, 1, base_length, vtkfile);

    /* seek back, write header block, seek forward */
    final_pos = ftell (vtkfile);
    base64_init_encodestate (&encode_state);
    base_length = base64_encode_block ((char *) compression_header, header_size, base_data, &encode_state);
    base_length += base64_encode_blockend (base_data + base_length, &encode_state);
    assert (base_length < code_length);
    base_data[base_length] = '\0';
    fseek1 = fseek (vtkfile, header_pos, SEEK_SET);
    (void) fwrite (base_data, 1, base_length, vtkfile);
    fseek2 = fseek (vtkfile, final_pos, SEEK_SET);

    /* clean up and return */
    free (compression_header);
    free (comp_data);
    free (base_data);
    if (fseek1 != 0 || fseek2 != 0 || ferror (vtkfile)) {
        return -1;
    }

    return 0;
}

#endif // defined(_USE_ZLIB) && defined(_USE_LIBB64)

#endif // KRATOS_VTK_H_INCLUDED

