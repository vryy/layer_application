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
#include "custom_utilities/spatial_grid_nodal_binning.h"
#include "custom_utilities/mesh_query_tool.h"
#include "custom_utilities/spatial_grid_elemental_binning.h"
#include "custom_utilities/structured_grid_elemental_indexing.h"
#include "custom_utilities/nonuniform_inclined_structured_grid_elemental_indexing.h"
#include "custom_utilities/model_part_utilities.h"
#include "custom_utilities/model_state.h"
#include "custom_python/add_utilities_to_python.h"

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

using namespace boost::python;

void Layer_AddTable(Layer& dummy, std::string table_name, boost::python::list& pyListNodes)
{
    std::vector<std::size_t> node_ids;
    typedef boost::python::stl_input_iterator<int> iterator_type;
    BOOST_FOREACH(const iterator_type::value_type& id,
                  std::make_pair(iterator_type(pyListNodes), // begin
                    iterator_type() ) ) // end
    {
        node_ids.push_back(id);
    }
    dummy.AddTable(table_name, node_ids);
}

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

template<class TEntityType>
void ModelPartUtilities_CalculateLocalSystem(ModelPartUtilities& rDummy,
    TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int& echo_level)
{
    rDummy.CalculateLocalSystem(rEntity, rCurrentProcessInfo, echo_level);
}

template<class TEntityType>
void ModelPartUtilities_CalculateMassMatrix(ModelPartUtilities& rDummy,
    TEntityType& rEntity, const ProcessInfo& rCurrentProcessInfo, const int& echo_level)
{
    rDummy.CalculateMassMatrix(rEntity, rCurrentProcessInfo, echo_level);
}

boost::python::list SpatialGridNodalBinning_GetNeighboursList(SpatialGridNodalBinning& dummy, ModelPart& r_model_part, int id, double r)
{
    boost::python::list list;
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
boost::python::list ElementalDataState_getitem(TElementalDataState& dummy, const std::size_t& key)
{
    typedef typename TElementalDataState::DataType DataType;
    const std::vector<DataType>& data = dummy[key];

    boost::python::list output;
    for (std::size_t i = 0; i < data.size(); ++i)
        output.append(data[i]);

    return output;
}

template<typename TVariableType>
void KratosLayerApplication_AddNodalDataStateToPython(const std::string& Name)
{
    typedef NodalDataState<TVariableType> NodalDataType;
    class_<NodalDataType, typename NodalDataType::Pointer, boost::noncopyable >
    ( Name.c_str(), init<TVariableType>() )
    .def("__getitem__", &NodalDataState_getitem<NodalDataType>)
    ;
}

template<typename TVariableType>
void KratosLayerApplication_AddElementalDataStateToPython(const std::string& Name)
{
    typedef ElementalDataState<TVariableType> ElementalDataType;
    class_<ElementalDataType, typename ElementalDataType::Pointer, boost::noncopyable >
    ( Name.c_str(), init<TVariableType>() )
    .def("__getitem__", &ElementalDataState_getitem<ElementalDataType>)
    ;
}

#ifdef LAYER_APP_USE_MMG
template<int TDim>
void LayerApp_ExportMMGMesher()
{
    void(MMGMesher<TDim>::*pointer_to_SetValueInt)(const Variable<int>&, const int&) = &MMGMesher<TDim>::SetValue;
    void(MMGMesher<TDim>::*pointer_to_SetParamInt)(const int&, const int&) = &MMGMesher<TDim>::SetValue;
    void(MMGMesher<TDim>::*pointer_to_SetValueDouble)(const Variable<double>&, const double&) = &MMGMesher<TDim>::SetValue;
    void(MMGMesher<TDim>::*pointer_to_SetParamDouble)(const int&, const double&) = &MMGMesher<TDim>::SetValue;

    std::stringstream ss;
    ss << "MMG" << TDim << "DMesher";
    class_<MMGMesher<TDim>, typename MMGMesher<TDim>::Pointer, boost::noncopyable>(ss.str().c_str(), init<>())
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

void LayerApp_AddCustomUtilitiesToPython()
{
    typedef ParameterList<std::string> ParameterListType;

    class_<Layer, Layer::Pointer, boost::noncopyable, bases<ParameterListType> >
    ("Layer", init<std::size_t, std::string>())
    .def("AddTable", &Layer_AddTable)
    .def(self_ns::str(self))
    ;

    class_<MDPAWriter, MDPAWriter::Pointer, boost::noncopyable>
    ("MDPAWriter", init<>())
    .def("WriteMDPA", &MDPAWriter::WriteMDPA)
    ;

    class_<MDPAModelPartWriter, MDPAModelPartWriter::Pointer, boost::noncopyable, bases<MDPAWriter> >
    ("MDPAModelPartWriter", init<const ModelPart&>())
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

    class_<AutoCollapseSpatialBinning<>, typename AutoCollapseSpatialBinning<>::Pointer, boost::noncopyable>
    ("AutoCollapseSpatialBinning", init<double, double, double, double, double, double, double>())
    .def("AddNode", &AutoCollapseSpatialBinning<>::AddNode)
    .def("NumberOfNodes", &AutoCollapseSpatialBinning<>::NumberOfNodes)
    .def("GetX", &AutoCollapseSpatialBinning<>::GetX)
    .def("GetY", &AutoCollapseSpatialBinning<>::GetY)
    .def("GetZ", &AutoCollapseSpatialBinning<>::GetZ)
    .def("SetDistance", &AutoCollapseSpatialBinning<>::SetDistance)
    .def("SetTolerance", &AutoCollapseSpatialBinning<>::SetTolerance)
    ;

    class_<SpatialGridNodalBinning, SpatialGridNodalBinning::Pointer, boost::noncopyable>
    ("SpatialGridNodalBinning", init<double, double, double, double, double, double, double>())
    .def("AddNodes", &SpatialGridNodalBinning::AddNodes)
    .def("GetNeighboursList", &SpatialGridNodalBinning_GetNeighboursList)
    ;

    class_<MeshQueryTool<Element>, typename MeshQueryTool<Element>::Pointer, boost::noncopyable>
    ("ElementalQueryTool", init<>())
    .def("SetTolerance", &MeshQueryTool<Element>::SetTolerance)
    .def("Initialize", &MeshQueryTool<Element>::Initialize)
    ;

    class_<SpatialGridElementalBinning, SpatialGridElementalBinning::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("SpatialGridElementalBinning", init<double, double, double, double, double, double, double>())
    ;

    class_<StructuredGridElementalIndexing<1>, typename StructuredGridElementalIndexing<1>::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("StructuredGridElementalIndexing1D", init<double, double, std::size_t>())
    ;

    class_<StructuredGridElementalIndexing<2>, typename StructuredGridElementalIndexing<2>::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("StructuredGridElementalIndexing2D", init<double, double, double, double, std::size_t, std::size_t>())
    ;

    class_<StructuredGridElementalIndexing<3>, typename StructuredGridElementalIndexing<3>::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("StructuredGridElementalIndexing3D", init<double, double, double, double, double, double, std::size_t, std::size_t, std::size_t>())
    ;

    class_<NonuniformInclinedStructuredGridElementalIndexing<1>, typename NonuniformInclinedStructuredGridElementalIndexing<1>::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("NonuniformInclinedStructuredGridElementalIndexing1D", init<const double>())
    .def("SetAxis", &NonuniformInclinedStructuredGridElementalIndexing<1>::SetAxis)
    ;

    class_<NonuniformInclinedStructuredGridElementalIndexing<2>, typename NonuniformInclinedStructuredGridElementalIndexing<2>::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("NonuniformInclinedStructuredGridElementalIndexing2D", init<const double>())
    .def("SetAxis", &NonuniformInclinedStructuredGridElementalIndexing<2>::SetAxis)
    ;

    class_<NonuniformInclinedStructuredGridElementalIndexing<3>, typename NonuniformInclinedStructuredGridElementalIndexing<3>::Pointer, bases<MeshQueryTool<Element> >, boost::noncopyable>
    ("NonuniformInclinedStructuredGridElementalIndexing3D", init<const double>())
    .def("SetAxis", &NonuniformInclinedStructuredGridElementalIndexing<3>::SetAxis)
    ;

    class_<ModelPartUtilities, ModelPartUtilities::Pointer, boost::noncopyable>
    ("ModelPartUtilities", init<>())
    .def("ExportNodesAndEdgesInfomation", &ModelPartUtilities_ExportNodesAndEdgesInfomation)
    .def("ExportNodesAndEdgesInfomationToGiD", &ModelPartUtilities_ExportNodesAndEdgesInfomationToGiD)
    .def("CreateElementFromCondition", &ModelPartUtilities_CreateElementFromCondition)
    .def("CreateElementFromNodes", &ModelPartUtilities_CreateElementFromNodes)
    .def("CalculateLocalSystem", &ModelPartUtilities_CalculateLocalSystem<Element>)
    .def("CalculateLocalSystem", &ModelPartUtilities_CalculateLocalSystem<Condition>)
    .def("CalculateMassMatrix", &ModelPartUtilities_CalculateMassMatrix<Element>)
    .def("CalculateMassMatrix", &ModelPartUtilities_CalculateMassMatrix<Condition>)
    ;

    //    KratosLayerApplication_AddNodalDataStateToPython<Variable<bool> >("BoolNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<int> >("IntegerNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<double> >("DoubleNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<array_1d<double, 3> > >("Array1DNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<Vector> >("VectorNodalDataState");
    KratosLayerApplication_AddNodalDataStateToPython<Variable<Matrix> >("MatrixNodalDataState");

//    KratosLayerApplication_AddElementalDataStateToPython<Variable<bool> >("BoolElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<int> >("IntegerElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<double> >("DoubleElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<array_1d<double, 3> > >("Array1DElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<Vector> >("VectorElementalDataState");
    KratosLayerApplication_AddElementalDataStateToPython<Variable<Matrix> >("MatrixElementalDataState");

    class_<ModelState, ModelState::Pointer, boost::noncopyable>
    ( "ModelState", init<>() )
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
    .def( "StoreNodalState", &ModelState::StoreNodalState<Variable<int> > )
    .def( "StoreNodalState", &ModelState::StoreNodalState<Variable<double> > )
    .def( "StoreNodalState", &ModelState::StoreNodalState<Variable<array_1d<double, 3> > > )
    .def( "StoreNodalState", &ModelState::StoreNodalState<Variable<Vector> > )
    .def( "StoreNodalState", &ModelState::StoreNodalState<Variable<Matrix> > )
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

    class_<HDF5PostUtility, HDF5PostUtility::Pointer, boost::noncopyable>("HDF5PostUtility", init<const std::string&>())
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
    enum_<MMG2D_Param>("MMG2D_Param")
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

    enum_<MMG3D_Param>("MMG3D_Param")
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

    LayerApp_ExportMMGMesher<2>();
    LayerApp_ExportMMGMesher<3>();
    #endif
}

}  // namespace Python.

} // Namespace Kratos
