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
#include "custom_io/vtk_mesh_container.h"
#include "custom_python/add_io_to_python.h"

namespace Kratos
{

namespace Python
{

typedef SDGidPostIO<GidSDIntegrationPointsContainer, GidSDMeshContainer> SDGidPostIOType;
typedef SDTikzPostIO<TikzIntegrationPointsContainer, TikzMeshContainer> SDTikzPostIOType;
typedef VtkIO<VtkMeshContainer> VtkIOType;

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

void (SDGidPostIOType::*pointer_to_flag_write_nodal_results)( const char* FlagName, Flags const& rFlag,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_bool_write_nodal_results)( Variable<bool> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_double_write_nodal_results)( Variable<double> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_array1d_write_nodal_results)( Variable<array_1d<double, 3> > const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_vector_write_nodal_results)( Variable<Vector> const& rVariable,
        SDGidPostIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &SDGidPostIOType::WriteNodalResults;
void (SDGidPostIOType::*pointer_to_matrix_write_nodal_results)( Variable<Matrix> const& rVariable,
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

void  LayerApplication_AddIOToPython()
{
    using namespace boost::python;

    class_<SDGidPostIOType, SDGidPostIOType::Pointer, boost::noncopyable>(
        "SDGidPostIO", init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>()
    )
    .def("WriteMesh", SDGidPostIO_WriteMesh)
    .def("WriteNodeMesh", SDGidPostIO_WriteNodeMesh)

    .def("InitializeMesh", &SDGidPostIOType::InitializeMesh)
    .def("FinalizeMesh", &SDGidPostIOType::FinalizeMesh)

    .def("InitializeResults", &SDGidPostIOType::InitializeResults)
    .def("FinalizeResults", &SDGidPostIOType::FinalizeResults)

    .def("WriteNodalResults", pointer_to_flag_write_nodal_results)
    .def("WriteNodalResults", pointer_to_bool_write_nodal_results)
    .def("WriteNodalResults", pointer_to_double_write_nodal_results)
    .def("WriteNodalResults", pointer_to_array1d_write_nodal_results)
    .def("WriteNodalResults", pointer_to_vector_write_nodal_results)
    .def("WriteNodalResults", pointer_to_matrix_write_nodal_results)

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
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

    class_<SDTikzPostIOType, SDTikzPostIOType::Pointer, boost::noncopyable>(
        "SDTikzPostIO", init<std::string const&, WriteDeformedMeshFlag, WriteConditionsFlag>()
    )
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

    class_<VtkIOType, VtkIOType::Pointer, boost::noncopyable>(
        "VtkIO", init<std::string const&, VTK_PostMode>()
    )
    .def("Initialize", &VtkIOType::Initialize)
    .def("Finalize", &VtkIOType::Finalize)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_double)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_array1d)
    .def("RegisterNodalResults", pointer_to_register_nodal_results_vector)
    ;
}

}  // namespace Python.

} // Namespace Kratos

