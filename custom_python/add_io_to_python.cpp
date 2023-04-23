/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

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

void IntegerPrintOnGaussPoints( SDGidPostIOType& dummy, const Variable<int>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
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

void  LayerApplication_AddIOToPython()
{
    using namespace boost::python;

    class_<SDGidPostIOType, SDGidPostIOType::Pointer, boost::noncopyable>
    ("SDGidPostIO", init<std::string const&, GiD_PostMode, MultiFileFlag, WriteDeformedMeshFlag, WriteConditionsFlag>())

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
    .def("PrintOnGaussPoints", IntegerPrintOnGaussPoints)
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

    enum_<VTK_PostMode>("VTKPostMode")
    .value("VTK_PostAscii", VTK_PostAscii)
    .value("VTK_PostBinary", VTK_PostBinary)
    ;

    enum_<VTK_PostFileFormat>("VTK_PostFileFormat")
    .value("VTK_PostVTU", VTK_PostVTU)
    .value("VTK_PostVTM", VTK_PostVTM)
    ;

    class_<VtkIOType, VtkIOType::Pointer, boost::noncopyable>
    ("VtkIO", init<std::string const&, const VTK_PostMode&>())
    .def("Initialize", &VtkIOType::Initialize)
    .def("Finalize", &VtkIOType::Finalize)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_double)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_array1d)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_vector)
    ;

    class_<VtkVTUIOType, VtkVTUIOType::Pointer, bases<VtkIOType>, boost::noncopyable>
    ("VtkVTUIO", init<std::string const&, const VTK_PostMode&>())
    ;

    class_<VtkVTMIOType, VtkVTMIOType::Pointer, bases<VtkIOType>, boost::noncopyable>
    ("VtkVTMIO", init<std::string const&, const VTK_PostMode&>())
    ;

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

