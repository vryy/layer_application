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
#include "custom_io/explicit_gid_io.h"
#include "custom_io/gid_integration_point_container.h"
#include "includes/gid_mesh_container.h"
#include "custom_python/add_io_to_python.h"

namespace Kratos
{

namespace Python
{

typedef ExplicitGidIO<GidIntegrationPointsContainer, GidMeshContainer> GidIOType;

void ExplicitGidIO_WriteNodeMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing nodes mesh "<<std::endl;
    dummy.WriteNodeMesh( rThisMesh );
    std::cout<<"end printing nodes mesh "<<std::endl;
}

void ExplicitGidIO_WriteMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing mesh "<<std::endl;
    dummy.WriteMesh( rThisMesh );
    std::cout<<"end printing mesh "<<std::endl;
}

void (GidIOType::*pointer_to_bool_write_nodal_results)( Variable<bool> const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIOType::WriteNodalResults;
void (GidIOType::*pointer_to_double_write_nodal_results)( Variable<double> const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIOType::WriteNodalResults;
void (GidIOType::*pointer_to_array1d_write_nodal_results)( Variable<array_1d<double, 3> > const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIOType::WriteNodalResults;
void (GidIOType::*pointer_to_vector_write_nodal_results)( Variable<Vector> const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIOType::WriteNodalResults;
void (GidIOType::*pointer_to_matrix_write_nodal_results)( Variable<Matrix> const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIOType::WriteNodalResults;

/////////////////////////////////////////////////////////////
/// NON-HISTORICAL DATABASE                               ///
////////////////////////////////////////////////////////////
void (GidIOType::*pointer_to_bool_write_nodal_results_NH)( Variable<bool> const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIOType::WriteNodalResultsNonHistorical;
void (GidIOType::*pointer_to_double_write_nodal_results_NH)( Variable<double> const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIOType::WriteNodalResultsNonHistorical;
void (GidIOType::*pointer_to_array1d_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIOType::WriteNodalResultsNonHistorical;
void (GidIOType::*pointer_to_matrix_write_nodal_results_NH)( Variable<Matrix > const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIOType::WriteNodalResultsNonHistorical;
void (GidIOType::*local_axes_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable,
        GidIOType::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIOType::WriteLocalAxesOnNodesNonHistorical;

void DoublePrintOnGaussPoints( GidIOType& dummy, const Variable<double>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void Array1DPrintOnGaussPoints( GidIOType& dummy, const Variable<array_1d<double,3> >& rVariable,
                                ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void VectorPrintOnGaussPoints( GidIOType& dummy, const Variable<Vector>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void MatrixPrintOnGaussPoints( GidIOType& dummy, const Variable<Matrix>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void  LayerApplication_AddIOToPython()
{
    using namespace boost::python;

    class_<GidIOType, GidIOType::Pointer, bases<IO>, boost::noncopyable>(
        "ExplicitGidIO",init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>())
    .def("WriteMesh",ExplicitGidIO_WriteMesh)
    .def("WriteNodeMesh",ExplicitGidIO_WriteNodeMesh)

    .def("InitializeMesh",&GidIOType::InitializeMesh)
    .def("FinalizeMesh",&GidIOType::FinalizeMesh)

    .def("InitializeResults",&GidIOType::InitializeResults)
    .def("FinalizeResults",&GidIOType::FinalizeResults)

    .def("WriteNodalResults",pointer_to_bool_write_nodal_results)
    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)
    .def("WriteNodalResults",pointer_to_vector_write_nodal_results)
    .def("WriteNodalResults",pointer_to_matrix_write_nodal_results)

//    .def("WriteLocalAxesOnNodes",local_axes_write_nodal_results)
    // NonHistorical
    .def("WriteNodalResultsNonHistorical",pointer_to_bool_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_double_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_array1d_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_matrix_write_nodal_results_NH)
    .def("WriteLocalAxesOnNodesNonHistorical",local_axes_write_nodal_results_NH)

    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)

    .def("Flush",&GidIOType::Flush)
    .def("ChangeOutputName",&GidIOType::ChangeOutputName)
    .def("CloseResultFile",&GidIOType::CloseResultFile)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

