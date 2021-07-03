/*
see layer_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: June 4, 2021 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_io/gid_binary_reader.h"
#include "custom_io/sd_gid_post_io.h"
#include "custom_io/sd_tikz_post_io.h"
#include "custom_io/gid_sd_integration_point_container.h"
#include "custom_io/gid_sd_mesh_container.h"
#include "custom_io/tikz_integration_point_container.h"
#include "custom_io/tikz_mesh_container.h"
#include "custom_io/vtk_io.h"
#include "custom_io/vtk_vtu_io.h"
#include "custom_io/vtk_vtm_io.h"
#include "custom_io/vtk_mesh_container.h"
#include "custom_python3/add_io_to_python.h"

namespace Kratos
{

namespace Python
{

typedef SDGidPostIO<GidSDIntegrationPointsContainer, GidSDMeshContainer> SDGidPostIOType;
typedef SDTikzPostIO<TikzIntegrationPointsContainer, TikzMeshContainer> SDTikzPostIOType;
typedef VtkIO<VtkMeshContainer> VtkIOType;
typedef VtkVTUIO<VtkMeshContainer> VtkVTUIOType;
typedef VtkVTMIO<VtkMeshContainer> VtkVTMIOType;

void SDGidPostIO_WriteNodeMesh( SDGidPostIOType& dummy, SDGidPostIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing nodes mesh "<<std::endl;
    dummy.WriteNodeMesh( rThisMesh );
    std::cout<<"end printing nodes mesh "<<std::endl;
}

void SDGidPostIO_WriteMesh( SDGidPostIOType& dummy, SDGidPostIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing mesh "<<std::endl;
    dummy.WriteMesh( rThisMesh );
    std::cout<<"end printing mesh "<<std::endl;
}

void (SDGidPostIOType::*pointer_to_flag_write_nodal_results1)( const char* FlagName, Flags const& rFlag,
        double SolutionTag ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_bool_write_nodal_results1)( Variable<bool> const& rVariable,
        double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_double_write_nodal_results1)( Variable<double> const& rVariable,
        double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_array1d_write_nodal_results1)( Variable<array_1d<double, 3> > const& rVariable,
        double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_vector_write_nodal_results1)( Variable<Vector> const& rVariable,
        double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_matrix_write_nodal_results1)( Variable<Matrix> const& rVariable,
        double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;

void (SDGidPostIOType::*pointer_to_flag_write_nodal_results2)( const char* FlagName, Flags const& rFlag,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_bool_write_nodal_results2)( Variable<bool> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_double_write_nodal_results2)( Variable<double> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_array1d_write_nodal_results2)( Variable<array_1d<double, 3> > const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_vector_write_nodal_results2)( Variable<Vector> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_matrix_write_nodal_results2)( Variable<Matrix> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;

/////////////////////////////////////////////////////////////
/// NON-HISTORICAL DATABASE                               ///
/////////////////////////////////////////////////////////////
void (SDGidPostIOType::*pointer_to_bool_write_nodal_results_NH)( Variable<bool> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
    = &SDGidPostIOType::WriteNodalResultsNonHistorical;
void (SDGidPostIOType::*pointer_to_double_write_nodal_results_NH)( Variable<double> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
    = &SDGidPostIOType::WriteNodalResultsNonHistorical;
void (SDGidPostIOType::*pointer_to_array1d_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
    = &SDGidPostIOType::WriteNodalResultsNonHistorical;
void (SDGidPostIOType::*pointer_to_matrix_write_nodal_results_NH)( Variable<Matrix > const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
    = &SDGidPostIOType::WriteNodalResultsNonHistorical;
void (SDGidPostIOType::*local_axes_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
    = &SDGidPostIOType::WriteLocalAxesOnNodesNonHistorical;

void PrintElementalPartitionIndex( SDGidPostIOType& dummy, const Variable<double>& rVariable,
                               ModelPart& r_model_part, double SolutionTag, int rank)
{
    dummy.PrintElementalPartitionIndex( rVariable, r_model_part, SolutionTag, 0, rank );
}

void DoublePrintOnGaussPoints( SDGidPostIOType& dummy, const Variable<double>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void Array1DPrintOnGaussPoints( SDGidPostIOType& dummy, const Variable<array_1d<double,3> >& rVariable,
                                ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void VectorPrintOnGaussPoints( SDGidPostIOType& dummy, const Variable<Vector>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void MatrixPrintOnGaussPoints( SDGidPostIOType& dummy, const Variable<Matrix>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

/////////////////////////////////////////////////////////////

void (VtkIOType::*pointer_to_register_nodal_results_double)( Variable<double> const&, const std::size_t& ) = &VtkIOType::RegisterNodalResults;
void (VtkIOType::*pointer_to_register_nodal_results_array1d)( Variable<array_1d<double, 3> > const&, const std::size_t& ) = &VtkIOType::RegisterNodalResults;
void (VtkIOType::*pointer_to_register_nodal_results_vector)( Variable<Vector> const&, const std::size_t&, const std::size_t& ) = &VtkIOType::RegisterNodalResults;

/////////////////////////////////////////////////////////////

pybind11::list GiDBinaryReader_GetMeshesName(GiDBinaryReader& rDummy)
{
    std::vector<std::string> mesh_list;
    rDummy.GetMeshesName(mesh_list);

    pybind11::list Output;
    for (std::size_t i = 0; i < mesh_list.size(); ++i)
        Output.append(mesh_list[i]);
    return Output;
}

pybind11::dict GiDBinaryReader_ReadNodalScalarValues(GiDBinaryReader& rDummy, const Variable<double>& rVariable, const std::size_t& step_index)
{
    std::map<std::size_t, double> Values;

    rDummy.ReadNodalScalarValues(rVariable.Name(), step_index, Values);

    pybind11::dict Output;
    for (auto it = Values.begin(); it != Values.end(); ++it)
    {
        std::string key = std::to_string(it->first);
        Output[key.c_str()] = it->second;
    }
    return Output;
}

pybind11::dict GiDBinaryReader_ReadNodalArray1DValues(GiDBinaryReader& rDummy, const Variable<array_1d<double, 3> >& rVariable, const std::size_t& step_index)
{
    std::map<std::size_t, std::vector<double> > Values;

    rDummy.ReadNodalVectorValues(rVariable.Name(), step_index, Values, 3);

    pybind11::dict Output;
    for (auto it = Values.begin(); it != Values.end(); ++it)
    {
        pybind11::list vector;
        for (std::size_t i = 0; i < it->second.size(); ++i)
            vector.append(it->second[i]);
        std::string key = std::to_string(it->first);
        Output[key.c_str()] = vector;
    }
    return Output;
}

pybind11::dict GiDBinaryReader_ReadNodalVectorValues(GiDBinaryReader& rDummy, const Variable<Vector>& rVariable, const std::size_t& step_index, const std::size_t& vector_size)
{
    std::map<std::size_t, std::vector<double> > Values;

    rDummy.ReadNodalVectorValues(rVariable.Name(), step_index, Values, vector_size);

    pybind11::dict Output;
    for (auto it = Values.begin(); it != Values.end(); ++it)
    {
        pybind11::list vector;
        for (std::size_t i = 0; i < it->second.size(); ++i)
            vector.append(it->second[i]);
        std::string key = std::to_string(it->first);
        Output[key.c_str()] = vector;
    }
    return Output;
}

/////////////////////////////////////////////////////////////

void  LayerApplication_AddIOToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<SDGidPostIOType, SDGidPostIOType::Pointer>
    (m, "SDGidPostIO")
    .def(init<std::string const&, GiD_PostMode, MultiFileFlag, WriteDeformedMeshFlag, WriteConditionsFlag>())

    .def("WriteMesh", SDGidPostIO_WriteMesh)
    .def("WriteNodeMesh", SDGidPostIO_WriteNodeMesh)

    .def("InitializeMesh", &SDGidPostIOType::InitializeMesh)
    .def("FinalizeMesh", &SDGidPostIOType::FinalizeMesh)

    .def("InitializeResults", &SDGidPostIOType::InitializeResults)
    .def("FinalizeResults", &SDGidPostIOType::FinalizeResults)

    .def("WriteNodalResults", pointer_to_flag_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_bool_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_double_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_array1d_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_vector_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_matrix_write_nodal_results1)

    .def("WriteNodalResults", pointer_to_flag_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_bool_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_double_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_array1d_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_vector_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_matrix_write_nodal_results2)

//    .def("WriteLocalAxesOnNodes",local_axes_write_nodal_results)
    // NonHistorical
    .def("WriteNodalResultsNonHistorical", pointer_to_bool_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical", pointer_to_double_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical", pointer_to_array1d_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical", pointer_to_matrix_write_nodal_results_NH)
    .def("WriteLocalAxesOnNodesNonHistorical", local_axes_write_nodal_results_NH)

    .def("PrintElementalPartitionIndex", PrintElementalPartitionIndex)
    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)

    .def("Flush", &SDGidPostIOType::Flush)
    .def("ChangeOutputName", &SDGidPostIOType::ChangeOutputName)
    .def("CloseResultFile", &SDGidPostIOType::CloseResultFile)
    .def("Reset", &SDGidPostIOType::Reset)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

    class_<SDTikzPostIOType, SDTikzPostIOType::Pointer>
    (m, "SDTikzPostIO")
    .def(init<std::string const&, WriteDeformedMeshFlag, WriteConditionsFlag>())
    .def("SetWriteId", &SDTikzPostIOType::SetWriteId)
    .def("SetCamera", &SDTikzPostIOType::SetCamera)
    .def("SetStyle", &SDTikzPostIOType::SetStyle)
    .def("SetCurrentStyle", &SDTikzPostIOType::SetCurrentStyle)
    .def("InitializeMesh", &SDTikzPostIOType::InitializeMesh)
    .def("FinalizeMesh", &SDTikzPostIOType::FinalizeMesh)
    .def("WriteMesh", &SDTikzPostIOType::WriteMesh)
    .def("WriteNodeMesh", &SDTikzPostIOType::WriteNodeMesh)
    ;

    enum_<VTK_PostMode>(m, "VTKPostMode")
    .value("VTK_PostAscii", VTK_PostAscii)
    .value("VTK_PostBinary", VTK_PostBinary)
    ;

    enum_<VTK_PostFileFormat>(m, "VTK_PostFileFormat")
    .value("VTK_PostVTU", VTK_PostVTU)
    .value("VTK_PostVTM", VTK_PostVTM)
    ;

    class_<VtkIOType, VtkIOType::Pointer>
    (m, "VtkIO")
    .def(init<std::string const&, const VTK_PostMode&>())
    .def("Initialize", &VtkIOType::Initialize)
    .def("Finalize", &VtkIOType::Finalize)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_double)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_array1d)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_vector)
    ;

    class_<VtkVTUIOType, VtkVTUIOType::Pointer, VtkIOType>
    (m, "VtkVTUIO")
    .def(init<std::string const&, const VTK_PostMode&>())
    ;

    class_<VtkVTMIOType, VtkVTMIOType::Pointer, VtkIOType>
    (m, "VtkVTMIO")
    .def(init<std::string const&, const VTK_PostMode&>())
    ;

    class_<GiDBinaryReader, GiDBinaryReader::Pointer>
    (m, "GiDBinaryReader")
    .def(init<std::string const&>())
    .def("GetMeshesName", &GiDBinaryReader_GetMeshesName)
    .def("ReadNodalScalarValues", &GiDBinaryReader_ReadNodalScalarValues)
    .def("ReadNodalArray1DValues", &GiDBinaryReader_ReadNodalArray1DValues)
    .def("ReadNodalVectorValues", &GiDBinaryReader_ReadNodalVectorValues)
    // .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

