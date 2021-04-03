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
#include "custom_utilities/selective_collapsible_layer_handler.h"
#include "custom_utilities/multipatch_layer_handler.h"
#include "custom_utilities/auto_collapse_spatial_binning.h"
#include "custom_utilities/spatial_grid_binning.h"
#include "custom_utilities/model_part_utilities.h"
#include "custom_python/add_utilities_to_python.h"

#ifdef LAYER_APP_USE_HDF5
#include "custom_utilities/hdf5_post_utility.h"
#endif

namespace Kratos
{

namespace Python
{

Layer::Pointer LayerHandler_getitem(LayerHandler& dummy, std::string name)
{
    return dummy[name];
}

void ModelPartUtilities_ExportNodesAndEdgesInfomation(ModelPartUtilities& rDummy,
    const std::string& filename, ModelPart& r_model_part,
    const bool& with_id, const bool& with_x, const bool& with_y, const bool& with_z,
    const std::string& separator, const int& precision)
{
    std::ofstream fid;
    fid.open(filename.c_str());

    rDummy.ExportNodalCoordinates(fid, r_model_part, with_id, with_x, with_y, with_z, separator, precision);
    rDummy.ExportEdgeInformation(fid, r_model_part.Elements(), separator);

    fid.close();
}

void ModelPartUtilities_ExportNodesAndEdgesInfomationToGiD(ModelPartUtilities& rDummy,
    const std::string& filename, ModelPart& r_model_part, const int& precision)
{
    std::ofstream fid;
    fid.open(filename.c_str());

    rDummy.ExportNodalCoordinatesToGiD(fid, r_model_part, precision);
    rDummy.ExportEdgeInformationToGiD(fid, r_model_part.Elements());

    fid.close();
}

Element::Pointer ModelPartUtilities_CreateElementFromCondition(ModelPartUtilities& rDummy,
    const std::string& sample_elem_name, const std::size_t& Id, Properties::Pointer pProperties,
    Condition::Pointer pCond)
{
    return rDummy.CreateEntity<Element>(sample_elem_name, Id, pProperties, pCond->GetGeometry());
}

Element::Pointer ModelPartUtilities_CreateElementFromNodes(ModelPartUtilities& rDummy,
    ModelPart& r_model_part, const std::string& sample_elem_name,
    const std::size_t& Id, Properties::Pointer pProperties, boost::python::list& node_ids)
{
    std::vector<std::size_t> node_list;
    typedef boost::python::stl_input_iterator<int> iterator_value_type;
    BOOST_FOREACH(const iterator_value_type::value_type& id,
                  std::make_pair(iterator_value_type(node_ids), // begin
                  iterator_value_type() ) ) // end
    {
        node_list.push_back(static_cast<std::size_t>(id));
    }

    return rDummy.CreateEntity<Element>(r_model_part, sample_elem_name, Id, pProperties, node_list);
}

void LayerApp_AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    typedef ParameterList<std::string> ParameterListType;

    class_<Layer, Layer::Pointer, boost::noncopyable, bases<ParameterListType> >
    ("Layer", init<std::size_t, std::string>())
    .def("AddTable", &Layer::AddTable)
    .def(self_ns::str(self))
    ;

    class_<MDPAWriter, MDPAWriter::Pointer, boost::noncopyable>
    ("MDPAWriter", init<>())
    .def("WriteMDPA", &MDPAWriter::WriteMDPA)
    ;

    class_<MDPAModelPartWriter, MDPAModelPartWriter::Pointer, boost::noncopyable, bases<MDPAWriter> >
    ("MDPAModelPartWriter", init<ModelPart::Pointer>())
    .def("SetNodeIndexOffset", &MDPAModelPartWriter::SetNodeIndexOffset)
    ;

    class_<LayerHandler, LayerHandler::Pointer, boost::noncopyable, bases<MDPAWriter> >
    ("LayerHandler", init<>())
    .def("__getitem__", &LayerHandler_getitem)
    .def("Has", &LayerHandler::Has)
    .def("AddLayer", &LayerHandler::AddLayer)
    .def("Check", &LayerHandler::Check)
    .def("RenumberAll", &LayerHandler::RenumberAll)
    .def("WriteLayers", &LayerHandler::WriteLayers)
    .def("WriteSimpleFunction", &LayerHandler::WriteSimpleFunction)
    .def(self_ns::str(self))
    ;

    class_<CollapsibleLayerHandler, CollapsibleLayerHandler::Pointer, boost::noncopyable, bases<LayerHandler> >
     ("CollapsibleLayerHandler", init<>())
    .def("SetSpacing", &CollapsibleLayerHandler::SetSpacing)
    .def("Collapse", &CollapsibleLayerHandler::Collapse)
    ;

    class_<SelectiveCollapsibleLayerHandler, SelectiveCollapsibleLayerHandler::Pointer, boost::noncopyable, bases<LayerHandler> >
     ("SelectiveCollapsibleLayerHandler", init<>())
    .def("SetSpacing", &SelectiveCollapsibleLayerHandler::SetSpacing)
    .def("AddGroup", &SelectiveCollapsibleLayerHandler::AddGroup)
    .def("CollapseLayer", &SelectiveCollapsibleLayerHandler::CollapseLayer)
    .def("CollapseGroup", &SelectiveCollapsibleLayerHandler::CollapseGroup)
    ;

    class_<MultipatchLayerHandler, MultipatchLayerHandler::Pointer, boost::noncopyable, bases<LayerHandler> >
     ("MultipatchLayerHandler", init<>())
    .def("AddLayerConnection", &MultipatchLayerHandler::AddLayerConnection)
    .def("FinalizeLayerConnection", &MultipatchLayerHandler::FinalizeLayerConnection)
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

    class_<ModelPartUtilities, ModelPartUtilities::Pointer, boost::noncopyable>
    ("ModelPartUtilities", init<>())
    .def("ExportNodesAndEdgesInfomation", &ModelPartUtilities_ExportNodesAndEdgesInfomation)
    .def("ExportNodesAndEdgesInfomationToGiD", &ModelPartUtilities_ExportNodesAndEdgesInfomationToGiD)
    .def("CreateElementFromCondition", &ModelPartUtilities_CreateElementFromCondition)
    .def("CreateElementFromNodes", &ModelPartUtilities_CreateElementFromNodes)
    ;

    #ifdef LAYER_APP_USE_HDF5
    class_<HDF5PostUtility, HDF5PostUtility::Pointer, boost::noncopyable>("HDF5PostUtility", init<const std::string&>())
    .def(init<const std::string&, const std::string&>())
    .def("WriteNodes", &HDF5PostUtility::WriteNodes)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<double>)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<array_1d<double, 3> >)
    .def("WriteNodalResults", &HDF5PostUtility::WriteNodalResults<Vector>)
    .def("WriteElementalData", &HDF5PostUtility::WriteElementalData<bool>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<double>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<array_1d<double, 3> >)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<Vector>)
    .def("ReadElementalData", &HDF5PostUtility::ReadElementalData<bool>)
    ;
    #endif
}

}  // namespace Python.

} // Namespace Kratos
