//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: 22 Oct 2015 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "custom_io/gid_binary_reader.h"
#include "custom_io/sd_gid_post_io.h"
#include "custom_io/sd_tikz_post_io.h"
#include "custom_io/gid_sd_integration_point_container.h"
#include "custom_io/gid_sd_mesh_container.h"
#include "custom_io/gid_mfem_mesh_container.h"
#include "custom_io/tikz_integration_point_container.h"
#include "custom_io/tikz_mesh_container.h"
#include "custom_io/vtk_io.h"
#include "custom_io/vtk_vtu_io.h"
#include "custom_io/vtk_vtm_io.h"
#include "custom_io/vtk_mesh_container.h"
#include "custom_python/add_io_to_python.h"

namespace Kratos
{

namespace Python
{

template<typename TGidPostIOType>
void GidPostIO_WriteNodeMesh( TGidPostIOType& dummy, typename TGidPostIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing nodes mesh "<<std::endl;
    dummy.WriteNodeMesh( rThisMesh );
    std::cout<<"end printing nodes mesh "<<std::endl;
}

template<typename TGidPostIOType>
void GidPostIO_WriteMesh( TGidPostIOType& dummy, typename TGidPostIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing mesh "<<std::endl;
    dummy.WriteMesh( rThisMesh );
    std::cout<<"end printing mesh "<<std::endl;
}

template<typename TGidPostIOType>
void PrintElementalPartitionIndex( TGidPostIOType& dummy, const Variable<double>& rVariable,
                               typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag, int rank)
{
    dummy.PrintElementalPartitionIndex( rVariable, r_model_part, SolutionTag, 0, rank );
}

template<typename TGidPostIOType>
void IntegerPrintOnGaussPoints( TGidPostIOType& dummy, const Variable<int>& rVariable,
                               typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

template<typename TGidPostIOType>
void DoublePrintOnGaussPoints( TGidPostIOType& dummy, const Variable<typename TGidPostIOType::DataType>& rVariable,
                               typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

template<typename TGidPostIOType>
void Array1DPrintOnGaussPoints( TGidPostIOType& dummy, const Variable<array_1d<typename TGidPostIOType::DataType, 3> >& rVariable,
                                typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

template<typename TGidPostIOType>
void VectorPrintOnGaussPoints1( TGidPostIOType& dummy, const Variable<typename TGidPostIOType::VectorType>& rVariable,
                               typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

template<typename TGidPostIOType>
void VectorPrintOnGaussPoints2( TGidPostIOType& dummy, const Variable<typename TGidPostIOType::VectorType>& rVariable,
                               typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag, int value_index )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag, value_index );
}

template<typename TGidPostIOType>
void MatrixPrintOnGaussPoints( TGidPostIOType& dummy, const Variable<typename TGidPostIOType::MatrixType>& rVariable,
                               typename TGidPostIOType::ModelPartType& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

/////////////////////////////////////////////////////////////

boost::python::list GiDBinaryReader_GetMeshesName(GiDBinaryReader& rDummy)
{
    std::vector<std::string> mesh_list;
    rDummy.GetMeshesName(mesh_list);

    boost::python::list Output;
    for (std::size_t i = 0; i < mesh_list.size(); ++i)
        Output.append(mesh_list[i]);
    return Output;
}

boost::python::dict GiDBinaryReader_ReadNodalScalarValues(GiDBinaryReader& rDummy, const Variable<double>& rVariable, const std::size_t& step_index)
{
    std::map<std::size_t, double> Values;

    rDummy.ReadNodalScalarValues(rVariable.Name(), step_index, Values);

    boost::python::dict Output;
    for (auto it = Values.begin(); it != Values.end(); ++it)
        Output[it->first] = it->second;
    return Output;
}

boost::python::dict GiDBinaryReader_ReadNodalArray1DValues(GiDBinaryReader& rDummy, const Variable<array_1d<double, 3> >& rVariable, const std::size_t& step_index)
{
    std::map<std::size_t, std::vector<double> > Values;

    rDummy.ReadNodalVectorValues(rVariable.Name(), step_index, Values, 3);

    boost::python::dict Output;
    for (auto it = Values.begin(); it != Values.end(); ++it)
    {
        boost::python::list vector;
        for (std::size_t i = 0; i < it->second.size(); ++i)
            vector.append(it->second[i]);
        Output[it->first] = vector;
    }
    return Output;
}

boost::python::dict GiDBinaryReader_ReadNodalVectorValues(GiDBinaryReader& rDummy, const Variable<Vector>& rVariable, const std::size_t& step_index, const std::size_t& vector_size)
{
    std::map<std::size_t, std::vector<double> > Values;

    rDummy.ReadNodalVectorValues(rVariable.Name(), step_index, Values, vector_size);

    boost::python::dict Output;
    for (auto it = Values.begin(); it != Values.end(); ++it)
    {
        boost::python::list vector;
        for (std::size_t i = 0; i < it->second.size(); ++i)
            vector.append(it->second[i]);
        Output[it->first] = vector;
    }
    return Output;
}

/////////////////////////////////////////////////////////////

template<typename TGidPostIOType>
void LayerApplication_AddGidPostIOToPython(const std::string& name)
{
    using namespace boost::python;

    void (TGidPostIOType::*pointer_to_flag_write_nodal_results1)( const char* FlagName, Flags const& rFlag,
            double SolutionTag ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_bool_write_nodal_results1)( Variable<bool> const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_int_write_nodal_results1)( Variable<int> const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_double_write_nodal_results1)( Variable<double> const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_complex_write_nodal_results1)( Variable<std::complex<double> > const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_array1d_write_nodal_results1)( Variable<array_1d<double, 3> > const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_complexarray1d_write_nodal_results1)( Variable<array_1d<std::complex<double>, 3> > const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_vector_write_nodal_results1)( Variable<Vector> const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_matrix_write_nodal_results1)( Variable<Matrix> const& rVariable,
            double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;

    void (TGidPostIOType::*pointer_to_flag_write_nodal_results2)( const char* FlagName, Flags const& rFlag,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_bool_write_nodal_results2)( Variable<bool> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_int_write_nodal_results2)( Variable<int> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_double_write_nodal_results2)( Variable<double> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_complex_write_nodal_results2)( Variable<std::complex<double> > const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_array1d_write_nodal_results2)( Variable<array_1d<double, 3> > const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_complexarray1d_write_nodal_results2)( Variable<array_1d<std::complex<double>, 3> > const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_vector_write_nodal_results2)( Variable<Vector> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;
    void (TGidPostIOType::*pointer_to_matrix_write_nodal_results2)( Variable<Matrix> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
            std::size_t SolutionStepNumber ) = &TGidPostIOType::WriteNodalResults;

    /////////////////////////////////////////////////////////////
    /// NON-HISTORICAL DATABASE                               ///
    /////////////////////////////////////////////////////////////
    void (TGidPostIOType::*pointer_to_bool_write_nodal_results_NH)( Variable<bool> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
        = &TGidPostIOType::WriteNodalResultsNonHistorical;
    void (TGidPostIOType::*pointer_to_double_write_nodal_results_NH)( Variable<double> const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
        = &TGidPostIOType::WriteNodalResultsNonHistorical;
    void (TGidPostIOType::*pointer_to_array1d_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
        = &TGidPostIOType::WriteNodalResultsNonHistorical;
    void (TGidPostIOType::*pointer_to_matrix_write_nodal_results_NH)( Variable<Matrix > const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
        = &TGidPostIOType::WriteNodalResultsNonHistorical;
    void (TGidPostIOType::*local_axes_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable,
            const typename TGidPostIOType::NodesContainerType& rNodes, double SolutionTag)
        = &TGidPostIOType::WriteLocalAxesOnNodesNonHistorical;

    class_<TGidPostIOType, typename TGidPostIOType::Pointer, boost::noncopyable>
    (name.c_str(), init<std::string const&, GiD_PostMode, MultiFileFlag, WriteDeformedMeshFlag, WriteConditionsFlag>())

    .def("WriteMesh", GidPostIO_WriteMesh<TGidPostIOType>)
    .def("WriteNodeMesh", GidPostIO_WriteNodeMesh<TGidPostIOType>)

    .def("InitializeMesh", &TGidPostIOType::InitializeMesh)
    .def("FinalizeMesh", &TGidPostIOType::FinalizeMesh)

    .def("InitializeResults", &TGidPostIOType::InitializeResults)
    .def("FinalizeResults", &TGidPostIOType::FinalizeResults)

    .def("WriteNodalResults", pointer_to_flag_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_bool_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_int_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_double_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_complex_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_array1d_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_complexarray1d_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_vector_write_nodal_results1)
    .def("WriteNodalResults", pointer_to_matrix_write_nodal_results1)

    .def("WriteNodalResults", pointer_to_flag_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_bool_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_int_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_double_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_complex_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_array1d_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_complexarray1d_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_vector_write_nodal_results2)
    .def("WriteNodalResults", pointer_to_matrix_write_nodal_results2)

//    .def("WriteLocalAxesOnNodes",local_axes_write_nodal_results)
    // NonHistorical
    .def("WriteNodalResultsNonHistorical", pointer_to_bool_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical", pointer_to_double_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical", pointer_to_array1d_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical", pointer_to_matrix_write_nodal_results_NH)
    .def("WriteLocalAxesOnNodesNonHistorical", local_axes_write_nodal_results_NH)

    .def("PrintElementalPartitionIndex", PrintElementalPartitionIndex<TGidPostIOType>)
    .def("PrintOnGaussPoints", IntegerPrintOnGaussPoints<TGidPostIOType>)
    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints<TGidPostIOType>)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPoints<TGidPostIOType>)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints1<TGidPostIOType>)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints2<TGidPostIOType>)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints<TGidPostIOType>)

    .def("Flush", &TGidPostIOType::Flush)
    .def("ChangeOutputName", &TGidPostIOType::ChangeOutputName)
    .def("CloseResultFile", &TGidPostIOType::CloseResultFile)
    .def("Reset", &TGidPostIOType::Reset)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;
}

void LayerApplication_AddIOToPython()
{
    using namespace boost::python;

    typedef SDGidPostIO<GidSDIntegrationPointsContainer<ModelPart>, GidSDMeshContainer<ModelPart> > SDGidPostIOType;
    typedef SDGidPostIO<GidSDIntegrationPointsContainer<ComplexModelPart>, GidSDMeshContainer<ComplexModelPart> > ComplexSDGidPostIOType;
    typedef SDGidPostIO<GidSDIntegrationPointsContainer<GComplexModelPart>, GidSDMeshContainer<GComplexModelPart, 1> > GComplexSDGidPostIORealType;
    typedef SDGidPostIO<GidSDIntegrationPointsContainer<GComplexModelPart>, GidSDMeshContainer<GComplexModelPart, 2> > GComplexSDGidPostIOImagType;
    typedef SDGidPostIO<GidSDIntegrationPointsContainer<ModelPart>, GidMfemMeshContainer<ModelPart> > MfemGidPostIOType;
    typedef SDTikzPostIO<TikzIntegrationPointsContainer, TikzMeshContainer> SDTikzPostIOType;
    typedef VtkIO<VtkMeshContainer> VtkIOType;
    typedef VtkVTUIO<VtkMeshContainer> VtkVTUIOType;
    typedef VtkVTMIO<VtkMeshContainer> VtkVTMIOType;

    /////////////////////////////////////////////////////////////

    LayerApplication_AddGidPostIOToPython<SDGidPostIOType>("SDGidPostIO");
    LayerApplication_AddGidPostIOToPython<ComplexSDGidPostIOType>("ComplexSDGidPostIO");
    LayerApplication_AddGidPostIOToPython<GComplexSDGidPostIORealType>("GComplexSDGidPostIOReal");
    LayerApplication_AddGidPostIOToPython<GComplexSDGidPostIOImagType>("GComplexSDGidPostIOImag");
    LayerApplication_AddGidPostIOToPython<MfemGidPostIOType>("MfemGidPostIO");

    /////////////////////////////////////////////////////////////

    class_<SDTikzPostIOType, SDTikzPostIOType::Pointer, boost::noncopyable>
    ("SDTikzPostIO", init<std::string const&, WriteDeformedMeshFlag, WriteConditionsFlag>())
    .def("SetWriteId", &SDTikzPostIOType::SetWriteId)
    .def("SetCamera", &SDTikzPostIOType::SetCamera)
    .def("SetStyle", &SDTikzPostIOType::SetStyle)
    .def("SetCurrentStyle", &SDTikzPostIOType::SetCurrentStyle)
    .def("InitializeMesh", &SDTikzPostIOType::InitializeMesh)
    .def("FinalizeMesh", &SDTikzPostIOType::FinalizeMesh)
    .def("WriteMesh", &SDTikzPostIOType::WriteMesh)
    .def("WriteNodeMesh", &SDTikzPostIOType::WriteNodeMesh)
    ;

    /////////////////////////////////////////////////////////////

    enum_<VTK_PostMode>("VTKPostMode")
    .value("VTK_PostAscii", VTK_PostAscii)
    .value("VTK_PostBinary", VTK_PostBinary)
    ;

    enum_<VTK_PostFileFormat>("VTK_PostFileFormat")
    .value("VTK_PostVTU", VTK_PostVTU)
    .value("VTK_PostVTM", VTK_PostVTM)
    ;

    void (VtkIOType::*pointer_to_register_nodal_results_double)( Variable<double> const&, const std::size_t ) = &VtkIOType::RegisterNodalResults;
    void (VtkIOType::*pointer_to_register_nodal_results_array1d)( Variable<array_1d<double, 3> > const&, const std::size_t ) = &VtkIOType::RegisterNodalResults;
    void (VtkIOType::*pointer_to_register_nodal_results_vector)( Variable<Vector> const&, const std::size_t, const std::size_t ) = &VtkIOType::RegisterNodalResults;
    void (VtkIOType::*pointer_to_register_cell_results_integer)( Variable<int> const&, const std::size_t ) = &VtkIOType::RegisterCellResults;
    void (VtkIOType::*pointer_to_register_cell_results_double)( Variable<double> const&, const std::size_t ) = &VtkIOType::RegisterCellResults;
    void (VtkIOType::*pointer_to_register_cell_properties)( ) = &VtkIOType::RegisterCellProperties;

    class_<VtkIOType, VtkIOType::Pointer, boost::noncopyable>
    ("VtkIO", init<std::string const&, const VTK_PostMode&>())
    .def("Initialize", &VtkIOType::Initialize)
    .def("Finalize", &VtkIOType::Finalize)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_double)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_array1d)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_vector)
    .def("RegisterCellResults", pointer_to_register_cell_results_integer)
    .def("RegisterCellResults", pointer_to_register_cell_results_double)
    .def("RegisterCellProperties", pointer_to_register_cell_properties)
    ;

    class_<VtkVTUIOType, VtkVTUIOType::Pointer, bases<VtkIOType>, boost::noncopyable>
    ("VtkVTUIO", init<std::string const&, const VTK_PostMode&>())
    ;

    class_<VtkVTMIOType, VtkVTMIOType::Pointer, bases<VtkIOType>, boost::noncopyable>
    ("VtkVTMIO", init<std::string const&, const VTK_PostMode&>())
    ;

    /////////////////////////////////////////////////////////////

    class_<GiDBinaryReader, GiDBinaryReader::Pointer, boost::noncopyable>
    ("GiDBinaryReader", init<std::string const&>())
    .def("GetMeshesName", &GiDBinaryReader_GetMeshesName)
    .def("ReadNodalScalarValues", &GiDBinaryReader_ReadNodalScalarValues)
    .def("ReadNodalArray1DValues", &GiDBinaryReader_ReadNodalArray1DValues)
    .def("ReadNodalVectorValues", &GiDBinaryReader_ReadNodalVectorValues)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

