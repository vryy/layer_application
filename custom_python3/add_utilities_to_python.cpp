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

// Project includes
#include "includes/define.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/mdpa_writer.h"
#include "custom_utilities/mdpa_model_part_writer.h"
#include "custom_utilities/layer.h"
// #include "custom_utilities/layer_handler.h"
// #include "custom_utilities/collapsible_layer_handler.h"
// #include "custom_utilities/selective_collapsible_layer_handler.h"
// #include "custom_utilities/multipatch_layer_handler.h"
#include "custom_utilities/auto_collapse_spatial_binning.h"
#include "custom_utilities/mesh_query_tool.h"
#include "custom_utilities/spatial_grid_elemental_binning.h"
#include "custom_utilities/spatial_grid_nodal_binning.h"
#include "custom_utilities/model_part_utilities.h"
#include "custom_utilities/model_state.h"
#include "custom_python3/add_utilities_to_python.h"

#ifdef LAYER_APP_USE_HDF5
#include "custom_utilities/hdf5_post_utility.h"
#endif

#ifdef LAYER_APP_USE_MMG
#include "custom_utilities/mmg_mesher.h"
#endif

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void Layer_AddTable(Layer& dummy, std::string table_name, pybind11::list& pyListNodes)
{
    std::vector<std::size_t> node_ids;
    for (std::size_t i = 0; i < pybind11::len(pyListNodes); ++i)
    {
        node_ids.push_back(static_cast<std::size_t>(pyListNodes[i].cast<int>()));
    }
    dummy.AddTable(table_name, node_ids);
}

// Layer::Pointer LayerHandler_getitem(LayerHandler& dummy, std::string name)
// {
//     return dummy[name];
// }

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
    const std::size_t& Id, Properties::Pointer pProperties, pybind11::list& node_ids)
{
    std::vector<std::size_t> node_list;
    for (std::size_t i = 0; i < pybind11::len(node_ids); ++i)
    {
        node_list.push_back(static_cast<std::size_t>(node_ids[i].cast<int>()));
    }

    return rDummy.CreateEntity<Element>(r_model_part, sample_elem_name, Id, pProperties, node_list);
}

template<class TEntityType>
void ModelPartUtilities_CalculateLocalSystem(ModelPartUtilities& rDummy,
    TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int& echo_level)
{
    rDummy.CalculateLocalSystem(rEntity, rCurrentProcessInfo, echo_level);
}

pybind11::list SpatialGridNodalBinning_GetNeighboursList(SpatialGridNodalBinning& dummy, ModelPart& r_model_part, int id, double r)
{
    pybind11::list list;
    auto Neighbours = dummy.GetNeighbourNodes(r_model_part, id, r);
    for(auto it = Neighbours.begin(); it != Neighbours.end(); ++it)
    {
        list.append(static_cast<int>(*it));
    }
    return list;
}

template<typename TNodalDataState>
typename TNodalDataState::DataType NodalDataState_getitem(TNodalDataState& dummy, const std::size_t& key)
{
    return dummy[key];
}

template<typename TElementalDataState>
pybind11::list ElementalDataState_getitem(TElementalDataState& dummy, const std::size_t& key)
{
    typedef typename TElementalDataState::DataType DataType;
    const std::vector<DataType>& data = dummy[key];

    pybind11::list output;
    for (std::size_t i = 0; i < data.size(); ++i)
        output.append(data[i]);

    return output;
}

template<typename TVariableType>
void KratosLayerApplication_AddNodalDataStateToPython(pybind11::module& m, const std::string& Name)
{
    typedef NodalDataState<TVariableType> NodalDataType;
    class_<NodalDataType, typename NodalDataType::Pointer>
    (m, Name.c_str())
    .def(init<const TVariableType&>())
    .def("__getitem__", &NodalDataState_getitem<NodalDataType>)
    ;
}

template<typename TVariableType>
void KratosLayerApplication_AddElementalDataStateToPython(pybind11::module& m, const std::string& Name)
{
    typedef ElementalDataState<TVariableType> ElementalDataType;
    class_<ElementalDataType, typename ElementalDataType::Pointer>
    (m, Name.c_str())
    .def(init<const TVariableType&>())
    .def("__getitem__", &ElementalDataState_getitem<ElementalDataType>)
    ;
}

#ifdef LAYER_APP_USE_MMG
template<int TDim>
void LayerApp_ExportMMGMesher(pybind11::module& m)
{
    void(MMGMesher<TDim>::*pointer_to_SetValueInt)(const Variable<int>&, const int&) = &MMGMesher<TDim>::SetValue;
    void(MMGMesher<TDim>::*pointer_to_SetParamInt)(const int&, const int&) = &MMGMesher<TDim>::SetValue;
    void(MMGMesher<TDim>::*pointer_to_SetValueDouble)(const Variable<double>&, const double&) = &MMGMesher<TDim>::SetValue;
    void(MMGMesher<TDim>::*pointer_to_SetParamDouble)(const int&, const double&) = &MMGMesher<TDim>::SetValue;

    std::stringstream ss;
    ss << "MMG" << TDim << "DMesher";
    class_<MMGMesher<TDim>, typename MMGMesher<TDim>::Pointer>
    (m, ss.str().c_str())
    .def(init<>())
    .def(init<const bool&, const bool&>())
    .def("SetValue", pointer_to_SetValueInt)
    .def("SetIParam", pointer_to_SetParamInt)
    .def("SetValue", pointer_to_SetValueDouble)
    .def("SetDParam", pointer_to_SetParamDouble)
    .def("Initialize", &MMGMesher<TDim>::Initialize)
    .def("Export", &MMGMesher<TDim>::Export)
    .def("SaveMesh", &MMGMesher<TDim>::SaveMesh)
    .def("SaveLevelSet", &MMGMesher<TDim>::SaveLevelSet)
    .def("SaveMetric", &MMGMesher<TDim>::SaveMetric)
    ;
}
#endif

void LayerApp_AddCustomUtilitiesToPython(pybind11::module& m)
{
    typedef ParameterList<std::string> ParameterListType;

    class_<Layer, Layer::Pointer, ParameterListType>
    (m, "Layer")
    .def(init<std::size_t, std::string>())
    .def("AddTable", &Layer_AddTable)
    // .def(self_ns::str(self))
    ;

    class_<MDPAWriter, MDPAWriter::Pointer>
    (m, "MDPAWriter")
    .def(init<>())
    .def("WriteMDPA", &MDPAWriter::WriteMDPA)
    ;

    class_<MDPAModelPartWriter, MDPAModelPartWriter::Pointer, MDPAWriter>
    (m, "MDPAModelPartWriter")
    .def(init<const ModelPart&>())
    .def("SetNodeIndexOffset", &MDPAModelPartWriter::SetNodeIndexOffset)
    ;

    // class_<LayerHandler, LayerHandler::Pointer, MDPAWriter>
    // ("LayerHandler", init<>())
    // .def("__getitem__", &LayerHandler_getitem)
    // .def("Has", &LayerHandler::Has)
    // .def("AddLayer", &LayerHandler::AddLayer)
    // .def("Check", &LayerHandler::Check)
    // .def("RenumberAll", &LayerHandler::RenumberAll)
    // .def("WriteLayers", &LayerHandler::WriteLayers)
    // .def("WriteSimpleFunction", &LayerHandler::WriteSimpleFunction)
    // .def(self_ns::str(self))
    // ;

    // class_<CollapsibleLayerHandler, CollapsibleLayerHandler::Pointer, LayerHandler>
    //  ("CollapsibleLayerHandler", init<>())
    // .def("SetSpacing", &CollapsibleLayerHandler::SetSpacing)
    // .def("Collapse", &CollapsibleLayerHandler::Collapse)
    // ;

    // class_<SelectiveCollapsibleLayerHandler, SelectiveCollapsibleLayerHandler::Pointer, LayerHandler>
    //  ("SelectiveCollapsibleLayerHandler", init<>())
    // .def("SetSpacing", &SelectiveCollapsibleLayerHandler::SetSpacing)
    // .def("AddGroup", &SelectiveCollapsibleLayerHandler::AddGroup)
    // .def("CollapseLayer", &SelectiveCollapsibleLayerHandler::CollapseLayer)
    // .def("CollapseGroup", &SelectiveCollapsibleLayerHandler::CollapseGroup)
    // ;

    // class_<MultipatchLayerHandler, MultipatchLayerHandler::Pointer, LayerHandler>
    //  ("MultipatchLayerHandler", init<>())
    // .def("AddLayerConnection", &MultipatchLayerHandler::AddLayerConnection)
    // .def("FinalizeLayerConnection", &MultipatchLayerHandler::FinalizeLayerConnection)
    // ;

    class_<AutoCollapseSpatialBinning, AutoCollapseSpatialBinning::Pointer>
    (m, "AutoCollapseSpatialBinning")
    .def(init<double, double, double, double, double, double, double>())
    .def("AddNode", &AutoCollapseSpatialBinning::AddNode)
    .def("NumberOfNodes", &AutoCollapseSpatialBinning::NumberOfNodes)
    .def("GetX", &AutoCollapseSpatialBinning::GetX)
    .def("GetY", &AutoCollapseSpatialBinning::GetY)
    .def("GetZ", &AutoCollapseSpatialBinning::GetZ)
    ;

    class_<SpatialGridNodalBinning, SpatialGridNodalBinning::Pointer>
    (m, "SpatialGridNodalBinning")
    .def(init<double, double, double, double, double, double, double>())
    .def("AddNodes", &SpatialGridNodalBinning::AddNodes)
    .def("GetNeighboursList", &SpatialGridNodalBinning_GetNeighboursList)
    ;

    class_<MeshQueryTool<Element>, typename MeshQueryTool<Element>::Pointer>
    (m, "ElementalQueryTool")
    .def(init<>())
    ;

    class_<SpatialGridElementalBinning, SpatialGridElementalBinning::Pointer>
    (m, "SpatialGridElementalBinning")
    .def(init<double, double, double, double, double, double, double>())
    ;

    class_<ModelPartUtilities, ModelPartUtilities::Pointer>
    (m, "ModelPartUtilities")
    .def(init<>())
    .def("ExportNodesAndEdgesInfomation", &ModelPartUtilities_ExportNodesAndEdgesInfomation)
    .def("ExportNodesAndEdgesInfomationToGiD", &ModelPartUtilities_ExportNodesAndEdgesInfomationToGiD)
    .def("CreateElementFromCondition", &ModelPartUtilities_CreateElementFromCondition)
    .def("CreateElementFromNodes", &ModelPartUtilities_CreateElementFromNodes)
    .def("CalculateLocalSystem", &ModelPartUtilities_CalculateLocalSystem<Element>)
    .def("CalculateLocalSystem", &ModelPartUtilities_CalculateLocalSystem<Condition>)
    ;

    //    KratosLayerApplication_AddNodalDataStateToPython<Variable<bool> >(m, "BoolNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<int> >(m, "IntegerNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<double> >(m, "DoubleNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<array_1d<double, 3> > >(m, "Array1DNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<Vector> >(m, "VectorNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<Matrix> >(m, "MatrixNodalDataState");

//    KratosLayerApplication_AddElementalDataStateToPython<Variable<bool> >(m, "BoolElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<int> >(m, "IntegerElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<double> >(m, "DoubleElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<array_1d<double, 3> > >(m, "Array1DElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<Vector> >(m, "VectorElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<Matrix> >(m, "MatrixElementalDataState");

    class_<ModelState, ModelState::Pointer>
    (m, "ModelState")
    .def(init<>())
//    .def( "GetNodalState", &ModelState::GetNodalState<Variable<bool> > )
    .def( "GetNodalState", &ModelState::GetNodalState<Variable<int> > )
    .def( "GetNodalState", &ModelState::GetNodalState<Variable<double> > )
    .def( "GetNodalState", &ModelState::GetNodalState<Variable<array_1d<double, 3> > > )
    .def( "GetNodalState", &ModelState::GetNodalState<Variable<Vector> > )
    .def( "GetNodalState", &ModelState::GetNodalState<Variable<Matrix> > )
//    .def( "GetElementalState", &ModelState::GetElementalState<Variable<bool> > )
    .def( "GetElementalState", &ModelState::GetElementalState<Variable<int> > )
    .def( "GetElementalState", &ModelState::GetElementalState<Variable<double> > )
    .def( "GetElementalState", &ModelState::GetElementalState<Variable<array_1d<double, 3> > > )
    .def( "GetElementalState", &ModelState::GetElementalState<Variable<Vector> > )
    .def( "GetElementalState", &ModelState::GetElementalState<Variable<Matrix> > )
//    .def( "SaveNodalState", &ModelState::SaveNodalState<Variable<bool> > )
    .def( "SaveNodalState", &ModelState::SaveNodalState<Variable<int> > )
    .def( "SaveNodalState", &ModelState::SaveNodalState<Variable<double> > )
    .def( "SaveNodalState", &ModelState::SaveNodalState<Variable<array_1d<double, 3> > > )
    .def( "SaveNodalState", &ModelState::SaveNodalState<Variable<Vector> > )
    .def( "SaveNodalState", &ModelState::SaveNodalState<Variable<Matrix> > )
//    .def( "SaveElementalState", &ModelState::SaveElementalState<Variable<bool> > )
    .def( "SaveElementalState", &ModelState::SaveElementalState<Variable<int> > )
    .def( "SaveElementalState", &ModelState::SaveElementalState<Variable<double> > )
    .def( "SaveElementalState", &ModelState::SaveElementalState<Variable<array_1d<double, 3> > > )
    .def( "SaveElementalState", &ModelState::SaveElementalState<Variable<Vector> > )
    .def( "SaveElementalState", &ModelState::SaveElementalState<Variable<Matrix> > )
//    .def( "LoadNodalState", &ModelState::LoadNodalState<Variable<bool> > )
    .def( "LoadNodalState", &ModelState::LoadNodalState<Variable<int> > )
    .def( "LoadNodalState", &ModelState::LoadNodalState<Variable<double> > )
    .def( "LoadNodalState", &ModelState::LoadNodalState<Variable<array_1d<double, 3> > > )
    .def( "LoadNodalState", &ModelState::LoadNodalState<Variable<Vector> > )
    .def( "LoadNodalState", &ModelState::LoadNodalState<Variable<Matrix> > )
//    .def( "LoadElementalState", &ModelState::LoadElementalState<Variable<bool> > )
    .def( "LoadElementalState", &ModelState::LoadElementalState<Variable<int> > )
    .def( "LoadElementalState", &ModelState::LoadElementalState<Variable<double> > )
    .def( "LoadElementalState", &ModelState::LoadElementalState<Variable<array_1d<double, 3> > > )
    .def( "LoadElementalState", &ModelState::LoadElementalState<Variable<Vector> > )
    .def( "LoadElementalState", &ModelState::LoadElementalState<Variable<Matrix> > )
    ;

    #ifdef LAYER_APP_USE_HDF5
    void(HDF5PostUtility::*pointer_to_WriteNodalResults_Nodes_double)(const Variable<double>&, const ModelPart::NodesContainerType&) = &HDF5PostUtility::WriteNodalResults<double>;
    void(HDF5PostUtility::*pointer_to_WriteNodalResults_Nodes_array_1d)(const Variable<array_1d<double, 3> >&, const ModelPart::NodesContainerType&) = &HDF5PostUtility::WriteNodalResults<array_1d<double, 3> >;
    void(HDF5PostUtility::*pointer_to_WriteNodalResults_Nodes_vector)(const Variable<Vector>&, const ModelPart::NodesContainerType&) = &HDF5PostUtility::WriteNodalResults<Vector>;
    void(HDF5PostUtility::*pointer_to_WriteNodalResults_ModelPart_double)(const Variable<double>&, ModelPart::Pointer) = &HDF5PostUtility::WriteNodalResults<double>;
    void(HDF5PostUtility::*pointer_to_WriteNodalResults_ModelPart_array_1d)(const Variable<array_1d<double, 3> >&, ModelPart::Pointer) = &HDF5PostUtility::WriteNodalResults<array_1d<double, 3> >;
    void(HDF5PostUtility::*pointer_to_WriteNodalResults_ModelPart_vector)(const Variable<Vector>&, ModelPart::Pointer) = &HDF5PostUtility::WriteNodalResults<Vector>;
    void(HDF5PostUtility::*pointer_to_WriteElementalData_ModelPart_bool)(const Variable<bool>&, ModelPart::Pointer) = &HDF5PostUtility::WriteElementalData<bool>;

    class_<HDF5PostUtility, HDF5PostUtility::Pointer>
    (m, "HDF5PostUtility")
    .def(init<const std::string&>())
    .def(init<const std::string&, const std::string&>())
    .def("WriteNodes", &HDF5PostUtility::WriteNodes)
    .def("WriteNodalResults", pointer_to_WriteNodalResults_Nodes_double)
    .def("WriteNodalResults", pointer_to_WriteNodalResults_Nodes_array_1d)
    .def("WriteNodalResults", pointer_to_WriteNodalResults_Nodes_vector)
    .def("WriteNodalResults", pointer_to_WriteNodalResults_ModelPart_double)
    .def("WriteNodalResults", pointer_to_WriteNodalResults_ModelPart_array_1d)
    .def("WriteNodalResults", pointer_to_WriteNodalResults_ModelPart_vector)
    .def("WriteElementalData", pointer_to_WriteElementalData_ModelPart_bool)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<double>)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<array_1d<double, 3> >)
    .def("ReadNodalResults", &HDF5PostUtility::ReadNodalResults<Vector>)
    .def("ReadElementalData", &HDF5PostUtility::ReadElementalData<bool>)
    ;
    #endif

    #ifdef LAYER_APP_USE_MMG
    enum_<MMG2D_Param>(m, "MMG2D_Param")
    .value("IPARAM_verbose", MMG2D_IPARAM_verbose)
    .value("IPARAM_mem", MMG2D_IPARAM_mem)
    .value("IPARAM_debug", MMG2D_IPARAM_debug)
    .value("IPARAM_angle", MMG2D_IPARAM_angle)
    .value("IPARAM_iso", MMG2D_IPARAM_iso)
    .value("IPARAM_opnbdy", MMG2D_IPARAM_opnbdy)
    .value("IPARAM_lag", MMG2D_IPARAM_lag)
    .value("IPARAM_3dMedit", MMG2D_IPARAM_3dMedit)
    .value("IPARAM_optim", MMG2D_IPARAM_optim)
    .value("IPARAM_noinsert", MMG2D_IPARAM_noinsert)
    .value("IPARAM_noswap", MMG2D_IPARAM_noswap)
    .value("IPARAM_nomove", MMG2D_IPARAM_nomove)
    .value("IPARAM_nosurf", MMG2D_IPARAM_nosurf)
    .value("IPARAM_nreg", MMG2D_IPARAM_nreg)
    .value("IPARAM_numsubdomain", MMG2D_IPARAM_numsubdomain)
    .value("IPARAM_numberOfLocalParam", MMG2D_IPARAM_numberOfLocalParam)
    .value("IPARAM_numberOfMat", MMG2D_IPARAM_numberOfMat)
    .value("IPARAM_anisosize", MMG2D_IPARAM_anisosize)
    .value("IPARAM_nosizreq", MMG2D_IPARAM_nosizreq)
    .value("DPARAM_angleDetection", MMG2D_DPARAM_angleDetection)
    .value("DPARAM_hmin", MMG2D_DPARAM_hmin)
    .value("DPARAM_hmax", MMG2D_DPARAM_hmax)
    .value("DPARAM_hsiz", MMG2D_DPARAM_hsiz)
    .value("DPARAM_hausd", MMG2D_DPARAM_hausd)
    .value("DPARAM_hgrad", MMG2D_DPARAM_hgrad)
    .value("DPARAM_hgradreq", MMG2D_DPARAM_hgradreq)
    .value("DPARAM_ls", MMG2D_DPARAM_ls)
    .value("DPARAM_rmc", MMG2D_DPARAM_rmc)
    ;

    enum_<MMG3D_Param>(m, "MMG3D_Param")
    .value("IPARAM_verbose", MMG3D_IPARAM_verbose)
    .value("IPARAM_mem", MMG3D_IPARAM_mem)
    .value("IPARAM_debug", MMG3D_IPARAM_debug)
    .value("IPARAM_angle", MMG3D_IPARAM_angle)
    .value("IPARAM_iso", MMG3D_IPARAM_iso)
    .value("IPARAM_nofem", MMG3D_IPARAM_nofem)
    .value("IPARAM_opnbdy", MMG3D_IPARAM_opnbdy)
    .value("IPARAM_lag", MMG3D_IPARAM_lag)
    .value("IPARAM_optim", MMG3D_IPARAM_optim)
    .value("IPARAM_optimLES", MMG3D_IPARAM_optimLES)
    .value("IPARAM_noinsert", MMG3D_IPARAM_noinsert)
    .value("IPARAM_noswap", MMG3D_IPARAM_noswap)
    .value("IPARAM_nomove", MMG3D_IPARAM_nomove)
    .value("IPARAM_nosurf", MMG3D_IPARAM_nosurf)
    .value("IPARAM_nreg", MMG3D_IPARAM_nreg)
    .value("IPARAM_numberOfLocalParam", MMG3D_IPARAM_numberOfLocalParam)
    .value("IPARAM_numberOfLSBaseReferences", MMG3D_IPARAM_numberOfLSBaseReferences)
    .value("IPARAM_numberOfMat", MMG3D_IPARAM_numberOfMat)
    .value("IPARAM_numsubdomain", MMG3D_IPARAM_numsubdomain)
    .value("IPARAM_renum", MMG3D_IPARAM_renum)
    .value("IPARAM_anisosize", MMG3D_IPARAM_anisosize)
    .value("IPARAM_octree", MMG3D_IPARAM_octree)
    .value("IPARAM_nosizreq", MMG3D_IPARAM_nosizreq)
    .value("DPARAM_angleDetection", MMG3D_DPARAM_angleDetection)
    .value("DPARAM_hmin", MMG3D_DPARAM_hmin)
    .value("DPARAM_hmax", MMG3D_DPARAM_hmax)
    .value("DPARAM_hsiz", MMG3D_DPARAM_hsiz)
    .value("DPARAM_hausd", MMG3D_DPARAM_hausd)
    .value("DPARAM_hgrad", MMG3D_DPARAM_hgrad)
    .value("DPARAM_hgradreq", MMG3D_DPARAM_hgradreq)
    .value("DPARAM_ls", MMG3D_DPARAM_ls)
    .value("DPARAM_rmc", MMG3D_DPARAM_rmc)
    .value("PARAM_size", MMG3D_PARAM_size)
    ;

    LayerApp_ExportMMGMesher<2>(m);
    LayerApp_ExportMMGMesher<3>(m);
    #endif
}

}  // namespace Python.

} // Namespace Kratos
