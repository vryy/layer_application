//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Oct 31, 2014 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/mdpa_writer.h"
#include "custom_utilities/mdpa_model_part_writer.h"
#include "custom_utilities/layer.h"
#include "custom_utilities/layer_handler.h"
#include "custom_utilities/collapsible_layer_handler.h"
#include "custom_utilities/multipatch_layer_handler.h"
#include "custom_utilities/auto_collapse_spatial_binning.h"
#include "custom_utilities/spatial_grid_binning.h"
#include "custom_python/add_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

Layer::Pointer LayerHandler_getitem(LayerHandler& dummy, std::string name)
{
    return dummy[name];
}

void LayerApp_AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    typedef ParameterList<std::string> ParameterListType;

    class_<Layer, Layer::Pointer, boost::noncopyable, bases<ParameterListType> >
    ("Layer", init<std::string>())
    .def("AddTable", &Layer::AddTable)
    .def(self_ns::str(self))
    ;

    class_<MDPAWriter, MDPAWriter::Pointer, boost::noncopyable>
    ("MDPAWriter", init<>())
    .def("WriteMDPA", &MDPAWriter::WriteMDPA)
    ;

    class_<MDPAModelPartWriter, MDPAModelPartWriter::Pointer, boost::noncopyable, bases<MDPAWriter> >
    ("MDPAModelPartWriter", init<ModelPart::Pointer>())
    ;

    class_<LayerHandler, LayerHandler::Pointer, boost::noncopyable, bases<MDPAWriter> >
    ("LayerHandler", init<>())
    .def("__getitem__", &LayerHandler_getitem)
    .def("Has", &LayerHandler::Has)
    .def("AddLayer", &LayerHandler::AddLayer)
    .def("RenumberAll", &LayerHandler::RenumberAll)
    .def("WriteLayers", &LayerHandler::WriteLayers)
    .def(self_ns::str(self))
    ;

    class_<CollapsibleLayerHandler, CollapsibleLayerHandler::Pointer, boost::noncopyable, bases<LayerHandler> >
     ("CollapsibleLayerHandler", init<>())
    .def("SetSpacing", &CollapsibleLayerHandler::SetSpacing)
    .def("AddGroup", &CollapsibleLayerHandler::AddGroup)
    .def("CollapseLayer", &CollapsibleLayerHandler::CollapseLayer)
    .def("CollapseGroup", &CollapsibleLayerHandler::CollapseGroup)
    ;

    class_<MultipatchLayerHandler, MultipatchLayerHandler::Pointer, boost::noncopyable, bases<LayerHandler> >
     ("MultipatchLayerHandler", init<>())
    .def("AddLayerConnection", &MultipatchLayerHandler::AddLayerConnection)
    ;

    class_<AutoCollapseSpatialBinning, AutoCollapseSpatialBinning::Pointer, boost::noncopyable>
    ("AutoCollapseSpatialBinning", init<double, double, double, double, double, double, double>())
    .def("AddNode", &AutoCollapseSpatialBinning::AddNode)
    .def("NumberOfNodes", &AutoCollapseSpatialBinning::NumberOfNodes)
    .def("GetX", &AutoCollapseSpatialBinning::GetX)
    .def("GetY", &AutoCollapseSpatialBinning::GetY)
    .def("GetZ", &AutoCollapseSpatialBinning::GetZ)
    ;

    class_<SpatialGridBinning, SpatialGridBinning::Pointer, boost::noncopyable>
    ("SpatialGridBinning", init<double, double, double, double, double, double, double>())
    .def("AddNodes", &SpatialGridBinning::AddNodes)
    .def("GetNeighboursList", &SpatialGridBinning::GetNeighboursList)
    ;

}

}  // namespace Python.

} // Namespace Kratos
