//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 31 Oct 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_LAYER_HANDLER_H_INCLUDED )
#define  KRATOS_LAYER_APP_LAYER_HANDLER_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <ctime>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "containers/pointer_vector_set.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/mdpa_writer.h"
#include "custom_utilities/layer.h"


namespace Kratos
{

///@addtogroup LayerApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/*** Detail class definition.
    Layer handler keeps a global node information and handle layer collapse properly
 */
class LayerHandler : public MDPAWriter
{
public:

    ///@name Type Definitions
    ///@{

    typedef Layer::IndexType IndexType;
    typedef Layer::PointType PointType;
    typedef Layer::EntityType EntityType;
    typedef PointerVectorSet<PointType, IndexedObject> PointsContainerType;
    typedef std::map<std::string, Layer::Pointer> LayersContainerType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(LayerHandler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LayerHandler()
    {
        mLastNode = 1;
        mLastEntity = 1;
        mLastPropId = 0;
        mLastLayer = 0;
    }

    /// Destructor.
    virtual ~LayerHandler()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    Layer::Pointer& operator[] (std::string name)
    {
        LayersContainerType::iterator it = mpLayers.find(name);
        if(it != mpLayers.end())
            return it->second;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Layer does not exist:", name)
    }

    Layer::Pointer& operator[] (std::size_t Id)
    {
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            if(it->second->Id() == Id)
            {
                return it->second;
            }
        }
        KRATOS_THROW_ERROR(std::logic_error, "Layer does not exist:", Id)
    }

    bool Has(std::string name) const
    {
        LayersContainerType::const_iterator it = mpLayers.find(name);
        return (it != mpLayers.end());
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void AddLayer(std::string layer_name, boost::python::dict& pyDictNodalSet,
            boost::python::dict& pyDictEntitySet, boost::python::dict& pyDictEntityInfoSet)
    {
        Layer::Pointer pLayer = Layer::Pointer(new Layer(++mLastLayer, layer_name));

        // iterate all the given nodes in layer and assign new nodal id; also build the map from the old node id to the new node id
        boost::python::list NodeIds = pyDictNodalSet.keys();
        typedef boost::python::stl_input_iterator<int> iterator_type;
        std::map<int, IndexType> NodeMap;
        BOOST_FOREACH(const iterator_type::value_type& id,
                      std::make_pair(iterator_type(NodeIds), // begin
                        iterator_type() ) ) // end
        {
            boost::python::object o = pyDictNodalSet.get(id);
            // here assumed that the nodal coordinates is given as a list
            boost::python::list coor = boost::python::extract<boost::python::list>(o);

            std::vector<double> P;
            typedef boost::python::stl_input_iterator<double> iterator_value_type;
            BOOST_FOREACH(const iterator_value_type::value_type& c,
                          std::make_pair(iterator_value_type(coor), // begin
                            iterator_value_type() ) ) // end
            {
                P.push_back(c);
            }

            PointType::Pointer pNewPoint = PointType::Pointer(new PointType(mLastNode, id, pLayer->Id(), P[0], P[1], P[2]));
            pLayer->AddNode(pNewPoint);
            mpPoints.push_back(pNewPoint);
            NodeMap[id] = mLastNode;

            ++mLastNode;
        }

        // iterate all the given entities in pLayer and reform the connectivities
        boost::python::list EntityIds = pyDictEntitySet.keys();
        std::map<int, IndexType> EntityMap;
        BOOST_FOREACH(const iterator_type::value_type& id,
                      std::make_pair(iterator_type(EntityIds), // begin
                        iterator_type() ) ) // end
        {
            boost::python::object o = pyDictEntitySet.get(id);
            // here assumed that the entity connectivities are given as a list
            boost::python::list connectivities = boost::python::extract<boost::python::list>(o);

            std::vector<IndexType> new_connectivities;
            EntityType::Pointer pNewEntity = EntityType::Pointer(new EntityType(mLastEntity));
            BOOST_FOREACH(const iterator_type::value_type& id,
                          std::make_pair(iterator_type(connectivities), // begin
                            iterator_type() ) ) // end
            {
                PointType::Pointer pPoint = mpPoints(NodeMap[id]);
                pNewEntity->AddNode(pPoint);
            }

            pLayer->AddEntity(pNewEntity);
            EntityMap.insert(std::pair<int, IndexType>(id, mLastEntity));
            ++mLastEntity;
        }

        // iterate all the given entity info set to set the value for each entity
        boost::python::list InfoNames = pyDictEntityInfoSet.keys();
        typedef boost::python::stl_input_iterator<std::string> iterator_info_type;
        BOOST_FOREACH(const iterator_info_type::value_type& info_name,
                      std::make_pair(iterator_info_type(InfoNames), // begin
                        iterator_info_type() ) ) // end
        {
            pLayer->AddInfo(info_name);
            // iterate all given entities associated with this info
            boost::python::object o = pyDictEntityInfoSet.get(info_name);
            boost::python::dict entities_info = boost::python::extract<boost::python::dict>(o);
            EntityIds = entities_info.keys();
            BOOST_FOREACH(const iterator_type::value_type& old_id,
                      std::make_pair(iterator_type(EntityIds), // begin
                        iterator_type() ) ) // end
            {
                IndexType new_id = EntityMap[old_id];
                boost::python::object o2 = entities_info.get(old_id);
                std::string o2_classname = boost::python::extract<std::string>(o2.attr("__class__").attr("__name__"));

                Layer::EntitiesIteratorType it = pLayer->find(new_id);
                if(it == pLayer->EntitiesEnd())
                    continue;
                if(o2_classname == std::string("int"))
                {
                    (*(*it))[info_name] = boost::python::extract<int>(o2);
                }
                else if(o2_classname == std::string("float"))
                {
                    (*(*it))[info_name] = boost::python::extract<double>(o2);
                }
                else if(o2_classname == std::string("list"))
                {
                    boost::python::list o3 = boost::python::extract<boost::python::list>(o2);
                    // extract the first element of o3 and check if it is a value or a list
                    boost::python::object o4 = boost::python::extract<boost::python::object>(o3[0]);
                    std::string o4_classname = boost::python::extract<std::string>(o4.attr("__class__").attr("__name__"));
                    if(o4_classname == std::string("float") || o4_classname == std::string("int"))
                    {
                        // a vector is detected
                        int len = boost::python::extract<int>(o3.attr("__len__")());
                        Vector v(len);
                        for(int i = 0; i < len; ++i)
                            v[i] = boost::python::extract<double>(o3[i]);
                        (*(*it))[info_name] = v;
                    }
                    else if(o4_classname == std::string("list"))
                    {
                        // a matrix is detected
                        int len1 = boost::python::extract<int>(o3.attr("__len__")());
                        int len2 = boost::python::extract<int>(o4.attr("__len__")());
                        Matrix M(len1, len2);
                        for(int i = 0; i < len1; ++i)
                        {
                            boost::python::list row = boost::python::extract<boost::python::list>(o3[i]);
                            int len_row = boost::python::extract<int>(row.attr("__len__")());
                            if(len_row != len2)
                                KRATOS_THROW_ERROR(std::logic_error, "Incompatiable matrix dimension detected at row", i)
                            else
                            {
                                for(int j = 0; j < len2; ++j)
                                    M(i, j) = boost::python::extract<double>(row[j]);
                            }
                        }
                        (*(*it))[info_name] = M;
                    }
                    else
                        KRATOS_THROW_ERROR(std::logic_error, "Unsupported underlying data type:", info_name)
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "Unsupported entity info type:", o2_classname)
            }
        }

        (*pLayer)["LAYER_PROP_ID"] = ++mLastPropId;
        mpLayers.insert(std::pair<std::string, Layer::Pointer>(layer_name, pLayer));
    }

    /// Check the integrity of all the layers
    void Check() const
    {
        for(LayersContainerType::const_iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);
            try
            {
                thisLayer.get<int>("LAYER_PROP_ID");
            }
            catch(std::exception& e)
            {
                std::stringstream ss;
                ss << "Layer " << thisLayer.Name() << " does not have a LAYER_PROP_ID";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        }
    }

    virtual void RenumberAll()
    {
        // collect all the unique id in the system
        std::set<IndexType> IdSet;
        for(PointsContainerType::iterator it = mpPoints.begin(); it != mpPoints.end(); ++it)
            IdSet.insert(it->Id());

        // for each unique id, assign the new id consecutively
        mLastNode = 0;
        std::map<IndexType, IndexType> NodeMap;
        for(std::set<IndexType>::iterator it = IdSet.begin(); it != IdSet.end(); ++it)
            NodeMap[*it] = ++mLastNode;

        // renumber all the nodes in nodes container
        for(PointsContainerType::iterator it = mpPoints.begin(); it != mpPoints.end(); ++it)
            it->SetId(NodeMap[it->Id()]);

        // renumber all the entities according to its type (element/condition)
        IndexType lastElement = 0;
        IndexType lastConditon = 0;
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

            for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
            {
                std::string layer_entity_type = thisLayer.get<std::string>("LAYER_ENTITY_TYPE");

                if(layer_entity_type == std::string("element") || layer_entity_type == std::string("bezier element")) // element
                {
                    (*it2)->SetId(++lastElement);
                }
                else if(layer_entity_type == std::string("condition") || layer_entity_type == std::string("bezier condition")) // condition
                {
                    (*it2)->SetId(++lastConditon);
                }
                else
                {
                    std::stringstream ss;
                    ss << "Unsupported entity type " << layer_entity_type << " of Layer " << thisLayer.Name();
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }
            }
        }

        std::cout << "RenumberAll completed" << std::endl;
    }

    void WriteLayers(std::string fn)
    {
        std::ofstream fid;
        std::string new_fn = fn + std::string(".py");
        fid.open(new_fn.c_str());

        // write element layers
        fid << "def ReadLayerSets():\n";
        fid << "\t" << "layer_sets = {}\n";
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

            std::string layer_entity_type = thisLayer.get<std::string>("LAYER_ENTITY_TYPE");

            if(layer_entity_type == std::string("element") || layer_entity_type == std::string("bezier element")) // element
            {}
            else if(layer_entity_type == std::string("condition") || layer_entity_type == std::string("bezier condition")) // condition
            {
                continue;
            }
            else if(layer_entity_type == std::string("null")) // none type
            {
                continue;
            }
            else
            {
                std::stringstream ss;
                ss << "Unsupported entity type " << layer_entity_type << " of Layer " << thisLayer.Name();
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            // start to write element layer
            fid << "\t" << "layer_sets['" << thisLayer.Name() << "'] = [";
            std::size_t cnt = 0;
            for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
            {
                fid << (*it2)->Id() << ", ";
                if(++cnt % 10 == 0)
                    fid << "\n\t";
            }
            fid << "]\n";
        }
        fid << "\t" << "return layer_sets\n\n";

        // write nodal layers
        fid << "def ReadLayerNodesSets():\n";
        fid << "\t" << "layer_nodes_sets = {}\n";
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

//            std::string layer_entity_type = thisLayer.get<std::string>("LAYER_ENTITY_TYPE");

//            if(layer_entity_type == std::string("element") || layer_entity_type == std::string("bezier element")) // element
//            {}
//            else if(layer_entity_type == std::string("condition") || layer_entity_type == std::string("bezier condition")) // condition
//            {}
//            else if(layer_entity_type == std::string("null")) // none type
//            {
//                continue;
//            }
//            else
//            {
//                std::stringstream ss;
//                ss << "Unsupported entity type " << layer_entity_type << " of Layer " << thisLayer.Name();
//                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
//            }

            // start to write node layer
            fid << "\t" << "layer_nodes_sets['" << thisLayer.Name() << "'] = [";
            std::size_t cnt = 0;
            for(Layer::NodesConstIteratorType it2 = thisLayer.NodesBegin(); it2 != thisLayer.NodesEnd(); ++it2)
            {
                fid << (*it2)->Id() << ", ";
                if(++cnt % 10 == 0)
                    fid << "\n\t";
            }
            fid << "]\n";
        }
        fid << "\t" << "return layer_nodes_sets\n\n";

        // write condition layers
        fid << "def ReadLayerConditionsSets():\n";
        fid << "\t" << "layer_conds_sets = {}\n";
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

            std::string layer_entity_type = thisLayer.get<std::string>("LAYER_ENTITY_TYPE");

            if(layer_entity_type == std::string("element") || layer_entity_type == std::string("bezier element")) // element
            {
                continue;
            }
            else if(layer_entity_type == std::string("condition") || layer_entity_type == std::string("bezier condition")) // condition
            {}
            else if(layer_entity_type == std::string("null")) // none type
            {
                continue;
            }
            else
            {
                std::stringstream ss;
                ss << "Unsupported entity type " << layer_entity_type << " of Layer " << thisLayer.Name();
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            // start to write condition layer
            fid << "\t" << "layer_conds_sets['" << thisLayer.Name() << "'] = [";
            std::size_t cnt = 0;
            for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
            {
                fid << (*it2)->Id() << ", ";
                if(++cnt % 10 == 0)
                    fid << "\n\t";
            }
            fid << "]\n";
        }
        fid << "\t" << "return layer_conds_sets\n\n";

        // write nodal tables
        fid << "def ReadTableNodesSets():\n";
        fid << "\t" << "table_nodes_sets = {}\n";
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);
            fid << "\t" << "table_nodes_sets['" << thisLayer.Name() << "'] = {}\n";

            std::vector<std::string> table_names = thisLayer.GetTableNames();
            for(std::size_t i = 0; i < table_names.size(); ++i)
            {
                fid << "\t" << "nodes_set = [";
                Layer::NodesContainerType table_nodes = thisLayer.Table(table_names[i]);
                std::size_t cnt = 0;
                for(Layer::NodesConstIteratorType it2 = table_nodes.ptr_begin(); it2 != table_nodes.ptr_end(); ++it2)
                {
                    fid << (*it2)->Id() << ", ";
                    if(++cnt % 10 == 0)
                        fid << "\n\t";
                }
                fid << "]\n";
                fid << "\t" << "table_nodes_sets['" << thisLayer.Name() << "']['" << table_names[i] << "'] = nodes_set\n";
            }
        }
        fid << "\t" << "return table_nodes_sets\n\n";

        // write the map from current node id to the original node id (ref id of node)
        fid << "def ReadLocalReferenceIds():\n";
        fid << "\t" << "layer_node_map = {}\n";
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

            fid << "\t" << "node_map = {}\n";
            for(Layer::NodesConstIteratorType it2 = thisLayer.NodesBegin(); it2 != thisLayer.NodesEnd(); ++it2)
            {
                fid << "\t" << "node_map[" << (*it2)->Id() << "] = " << (*it2)->LocalId() << "\n";
            }
            fid << "\t" << "layer_node_map['" << thisLayer.Name() << "'] = node_map\n";
        }
        fid << "\t" << "return layer_node_map\n\n";

        fid.close();
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "LayerHandler";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        std::cout << "+Nodes:" << std::endl;
        for(PointsContainerType::const_iterator it = mpPoints.begin(); it != mpPoints.end(); ++it)
            std::cout << "  " << it->Id() << ": " << it->X() << ", " << it->Y() << ", " << it->Z() << std::endl;
        for(LayersContainerType::const_iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
            std::cout << "+" << *(it->second);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    PointsContainerType mpPoints; // this container contains all the spatial points in the model
    LayersContainerType mpLayers; // layer container
    IndexType mLastNode;
    IndexType mLastEntity;
    int mLastPropId;
    IndexType mLastLayer;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void MDPA_Properties(std::ostream& rOStream)
    {
        typedef ParameterList<std::string> ParameterListType;
//        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
//        {
//            Layer& thisLayer = *(it->second);
//            int layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");
//            Begin(rOStream, "Properties", layer_prop_id);
//            for(ParameterListType::pair_iterator it2 = thisLayer.pair_begin(); it2 != thisLayer.pair_end(); ++it2)
//            {
//                // ignore the default underlying attributes
//                if(it2->first == std::string("LAYER_PROP_ID")
//                    || it2->first == std::string("LAYER_ENTITY_TYPE")
//                    || it2->first == std::string("LAYER_ENTITY_NAME"))
//                    continue;
//                rOStream << it2->first << " " << it2->second << std::endl;
//            }
//            End(rOStream, "Properties");
//        }

        // firstly collect all possible properties id
        std::set<int> prop_id_set;
        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);
            int layer_prop_id;
            try
            {
                layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");
            }
            catch(std::exception& e)
            {
                std::stringstream ss;
                ss << "Layer " << thisLayer.Name() << " does not have a LAYER_PROP_ID";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
            prop_id_set.insert(layer_prop_id);
        }

        for(std::set<int>::iterator it_id = prop_id_set.begin(); it_id != prop_id_set.end(); ++it_id)
        {
            Begin(rOStream, "Properties", *it_id);
            for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
            {
                Layer& thisLayer = *(it->second);
                int layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");
                if(layer_prop_id == *it_id)
                {
                    for(ParameterListType::pair_iterator it2 = thisLayer.pair_begin(); it2 != thisLayer.pair_end(); ++it2)
                    {
                        // ignore the default underlying attributes
                        if(it2->first == std::string("LAYER_PROP_ID")
                            || it2->first == std::string("LAYER_ENTITY_TYPE")
                            || it2->first == std::string("LAYER_ENTITY_NAME"))
                            continue;
                        rOStream << it2->first << " " << it2->second << std::endl;
                    }
                }
            }
            End(rOStream, "Properties");
        }
    }

    virtual void MDPA_Nodes(std::ostream& rOStream)
    {
        PointsContainerType pNodeList = mpPoints;
        pNodeList.Unique();

        Begin(rOStream, "Nodes");
        for(PointsContainerType::iterator it = pNodeList.begin(); it != pNodeList.end(); ++it)
        {
            rOStream << it->Id()
                     << " " << it->X()
                     << " " << it->Y()
                     << " " << it->Z()
                     << std::endl;
        }
        End(rOStream, "Nodes");
    }

    virtual void MDPA_Elements(std::ostream& rOStream)
    {
        typedef Layer::EntityType::BaseType::iterator EntityInfoIterator;

        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

            std::string layer_entity_type;
            try
            {
                layer_entity_type = thisLayer.get<std::string>("LAYER_ENTITY_TYPE");
            }
            catch(std::exception& e)
            {
                std::stringstream ss;
                ss << "Layer " << thisLayer.Name() << " does not have LAYER_ENTITY_TYPE";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            std::string layer_entity_name;
            try
            {
                layer_entity_name = thisLayer.get<std::string>("LAYER_ENTITY_NAME");
            }
            catch(std::exception& e)
            {
                std::stringstream ss;
                ss << "Layer " << thisLayer.Name() << " does not have LAYER_ENTITY_NAME";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            if(layer_entity_type == std::string("element")) // element
            {
                int layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");
//                KRATOS_WATCH(layer_entity_name)
//                KRATOS_WATCH(thisLayer.Name())

                // write the entity connectivities
                Begin(rOStream, "Elements", layer_entity_name);
                for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
                {
                    rOStream << (*it2)->Id() << " " << layer_prop_id;
                    for(std::size_t i = 0; i < (*it2)->size(); ++i)
                        rOStream << " " << (*it2)->operator[](i).Id();
                    rOStream << std::endl;
                }
                End(rOStream, "Elements");

                //write the entity info
                for(Layer::InfoNamesIteratorType it2 = thisLayer.InfoNamesBegin(); it2 != thisLayer.InfoNamesEnd(); ++it2)
                {
                    Begin(rOStream, "ElementalData", *it2);
                    for(Layer::EntitiesIteratorType it3 = thisLayer.EntitiesBegin(); it3 != thisLayer.EntitiesEnd(); ++it3)
                    {
                        EntityInfoIterator it4 = (*it3)->find(*it2);
                        // check if this entity contain the info
                        if(it4 != (*it3)->end())
                        {
                            rOStream << (*it3)->Id() << " " << *it4 << std::endl;
                        }
                    }
                    End(rOStream, "ElementalData");
                }
                std::cout << "Layer name " << thisLayer.Name() << " of type " << layer_entity_name << "(" << layer_entity_type << ") is written to Mdpa" << std::endl;
            }
            else if(layer_entity_type == std::string("bezier element")) // bezier element
            {
                int layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");

                // write the bezier block
                Begin(rOStream, "BezierBlock", "");

                rOStream << "    ";
                Begin(rOStream, "IsogeometricBezierData", "");

                std::vector<int> rowPtr;
                std::vector<int> colInd;
                std::vector<double> values;
                for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
                {
                    rOStream << "        "
                             << (*it2)->Id() << " " << (*it2)->size() << " ";

                    if(layer_entity_name.find("Bezier2D3") != std::string::npos)
                    {
                        rOStream << "2 3 ";
                    }
                    else if(layer_entity_name.find("Bezier2D") != std::string::npos)
                    {
                        rOStream << "2 2 ";
                    }
                    else if(layer_entity_name.find("Bezier3D") != std::string::npos)
                    {
                        rOStream << "3 3 ";
                    }
                    else if(layer_entity_name.find("Geo3dBezier") != std::string::npos)
                    {
                        rOStream << "3 3 ";
                    }
                    else
                    {
                        KRATOS_THROW_ERROR(std::runtime_error, "What's the hell Bezier", layer_entity_name)
                    }

                    Layer::InfoNamesIteratorType it3;
                    it3 = std::find(thisLayer.InfoNamesBegin(), thisLayer.InfoNamesEnd(), "NURBS_DEGREE_1");
                    int order_u = 0, order_v = 0, order_w = 0;
                    if(it3 != thisLayer.InfoNamesEnd())
                    {
                        EntityInfoIterator it4 = (*it2)->find(*it3);
                        // check if this entity contain the info
                        if(it4 != (*it2)->end())
                            order_u = boost::get<int>(*it4);
                    }
                    it3 = std::find(thisLayer.InfoNamesBegin(), thisLayer.InfoNamesEnd(), "NURBS_DEGREE_2");
                    if(it3 != thisLayer.InfoNamesEnd())
                    {
                        EntityInfoIterator it4 = (*it2)->find(*it3);
                        // check if this entity contain the info
                        if(it4 != (*it2)->end())
                            order_v = boost::get<int>(*it4);
                    }
                    it3 = std::find(thisLayer.InfoNamesBegin(), thisLayer.InfoNamesEnd(), "NURBS_DEGREE_3");
                    if(it3 != thisLayer.InfoNamesEnd())
                    {
                        EntityInfoIterator it4 = (*it2)->find(*it3);
                        // check if this entity contain the info
                        if(it4 != (*it2)->end())
                            order_w = boost::get<int>(*it4);
                    }

                    rOStream << order_u << " " << order_v << " " << order_w << std::endl;

                    //write the weights
                    EntityInfoIterator it4 = (*it2)->find("NURBS_WEIGHT");
                    if(it4 == (*it2)->end())
                        KRATOS_THROW_ERROR(std::runtime_error, "NURBS_WEIGHT is not provided for entity", (*it2)->Id())
                    rOStream << "        " << *it4 << std::endl;

                    //write the extraction operator

                    if((it4 = (*it2)->find("EXTRACTION_OPERATOR")) != (*it2)->end())
                    {
                        rOStream << "        Full" << std::endl;
                        rOStream << "        " << *it4 << std::endl;
                    }
                    else if((it4 = (*it2)->find("EXTRACTION_OPERATOR_MCSR")) != (*it2)->end())
                    {
                        rOStream << "        MCSR" << std::endl;
                        rOStream << "        " << *it4 << std::endl;
                    }
                    else if((it4 = (*it2)->find("EXTRACTION_OPERATOR_CSR")) != (*it2)->end())
                    {
//                        rOStream << "        CSR" << std::endl;
//                        rOStream << "        " << *it4 << std::endl;
                        KRATOS_THROW_ERROR(std::runtime_error, "EXTRACTION_OPERATOR_CSR is not implemented yet. Please wat.", "")
                    }
                    else
                        KRATOS_THROW_ERROR(std::runtime_error, "EXTRACTION_OPERATOR|EXTRACTION_OPERATOR_MCSR|EXTRACTION_OPERATOR_CSR is not provided for entity", (*it2)->Id())

                    rOStream << std::endl;
                }

                rOStream << "    ";
                End(rOStream, "IsogeometricBezierData");

                rOStream << "    ";
                Begin(rOStream, "ElementsWithGeometry", layer_entity_name);
                for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
                {
                    rOStream << "        " << (*it2)->Id() << " " << layer_prop_id << " " << (*it2)->Id();
                    for(unsigned int i = 0; i < (*it2)->size(); ++i)
                        rOStream << " " << (*it2)->operator[](i).Id();
                    rOStream << std::endl;
                }
                rOStream << "    ";
                End(rOStream, "ElementsWithGeometry");

                End(rOStream, "BezierBlock");

                //write the entity info
                for(Layer::InfoNamesIteratorType it2 = thisLayer.InfoNamesBegin(); it2 != thisLayer.InfoNamesEnd(); ++it2)
                {
                    if(it2->find("EXTRACTION_OPERATOR") != std::string::npos
                        || it2->find("NURBS_WEIGHT") != std::string::npos
                        || it2->find("NURBS_DEGREE") != std::string::npos )
                    {
                        continue;
                    }

                    Begin(rOStream, "ElementalData", *it2);
                    for(Layer::EntitiesIteratorType it3 = thisLayer.EntitiesBegin(); it3 != thisLayer.EntitiesEnd(); ++it3)
                    {
                        EntityInfoIterator it4 = (*it3)->find(*it2);
                        // check if this entity contain the info
                        if(it4 != (*it3)->end())
                        {
                            rOStream << (*it3)->Id() << " " << *it4 << std::endl;
                        }
                    }
                    End(rOStream, "ElementalData");
                }

                std::cout << "Layer name " << thisLayer.Name() << " of type " << layer_entity_name << "(" << layer_entity_type << ") is written to Mdpa" << std::endl;
            }
            else if(layer_entity_type == std::string("condition")) // condition
            {
                continue;
            }
            else if(layer_entity_type == std::string("bezier condition")) // condition
            {
                continue;
            }
            else if(layer_entity_type == std::string("null")) // null element/condition will not be accounted
            {
                continue;
            }
            else
            {
                std::stringstream ss;
                ss << "Unsupported entity type " << layer_entity_type << " of Layer " << thisLayer.Name();
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        }
    }

    virtual void MDPA_Conditions(std::ostream& rOStream)
    {
        typedef Layer::EntityType::BaseType::iterator EntityInfoIterator;

        for(LayersContainerType::iterator it = mpLayers.begin(); it != mpLayers.end(); ++it)
        {
            Layer& thisLayer = *(it->second);

            std::string layer_entity_type = thisLayer.get<std::string>("LAYER_ENTITY_TYPE");
//            KRATOS_WATCH(layer_entity_type)
            std::string layer_entity_name = thisLayer.get<std::string>("LAYER_ENTITY_NAME");
//            KRATOS_WATCH(layer_entity_name)

            std::string layer_entity_type_name;
            std::string layer_entity_info_type_name;
            if(layer_entity_type == std::string("condition")) // condition
            {
                layer_entity_type_name = std::string("Conditions");
                layer_entity_info_type_name = std::string("ConditionalData");

                int layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");

                // write the entity connectivities
                Begin(rOStream, layer_entity_type_name, layer_entity_name);
                for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
                {
                    rOStream << (*it2)->Id() << " " << layer_prop_id;
                    for(unsigned int i = 0; i < (*it2)->size(); ++i)
                        rOStream << " " << (*it2)->operator[](i).Id();
                    rOStream << std::endl;
                }
                End(rOStream, layer_entity_type_name);

                //write the entity info
                for(Layer::InfoNamesIteratorType it2 = thisLayer.InfoNamesBegin(); it2 != thisLayer.InfoNamesEnd(); ++it2)
                {
                    Begin(rOStream, layer_entity_info_type_name, *it2);
                    for(Layer::EntitiesIteratorType it3 = thisLayer.EntitiesBegin(); it3 != thisLayer.EntitiesEnd(); ++it3)
                    {
                        typedef Layer::EntityType::BaseType::iterator EntityInfoIterator;
                        EntityInfoIterator it4 = (*it3)->find(*it2);
                        // check if this entity contain the info
                        if(it4 != (*it3)->end())
                        {
                            rOStream << (*it3)->Id() << " " << *it4 << std::endl;
                        }
                    }
                    End(rOStream, layer_entity_info_type_name);
                }
                std::cout << "Layer name " << thisLayer.Name() << " of type " << layer_entity_name << "(" << layer_entity_type << ") is written to Mdpa" << std::endl;
            }
            else if(layer_entity_type == std::string("element")) // element
            {
                continue;
            }
            else if(layer_entity_type == std::string("bezier element")) // bezier element
            {
                continue;
            }
            else if(layer_entity_type == std::string("bezier condition")) // bezier condition
            {
                int layer_prop_id = thisLayer.get<int>("LAYER_PROP_ID");

                Begin(rOStream, "BezierBlock", "");

                rOStream << "    ";
                Begin(rOStream, "IsogeometricBezierData", "");

                std::vector<int> rowPtr;
                std::vector<int> colInd;
                std::vector<double> values;
                for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
                {
                    rOStream << "        "
                             << (*it2)->Id() << " " << (*it2)->size() << " ";

                    if(layer_entity_name.find("Bezier2D3") != std::string::npos
                        || layer_entity_name.find("FaceLoadBezier") != std::string::npos
                        || layer_entity_name.find("FacePressureBezier") != std::string::npos)
                    {
                        rOStream << "2 3 ";
                    }
                    else if(layer_entity_name.find("Bezier2D") != std::string::npos)
                    {
                        rOStream << "2 2 ";
                    }
                    else if(layer_entity_name.find("Bezier3D") != std::string::npos)
                    {
                        rOStream << "3 3 ";
                    }
                    else if(layer_entity_name.find("Geo3dBezier") != std::string::npos)
                    {
                        rOStream << "3 3 ";
                    }
                    else
                    {
                        KRATOS_THROW_ERROR(std::runtime_error, "What's the fuck Bezier", layer_entity_name)
                    }

                    Layer::InfoNamesIteratorType it3;
                    it3 = std::find(thisLayer.InfoNamesBegin(), thisLayer.InfoNamesEnd(), "NURBS_DEGREE_1");
                    int order_u = 0, order_v = 0, order_w = 0;
                    if(it3 != thisLayer.InfoNamesEnd())
                    {
                        EntityInfoIterator it4 = (*it2)->find(*it3);
                        // check if this entity contain the info
                        if(it4 != (*it2)->end())
                            order_u = boost::get<int>(*it4);
                    }
                    it3 = std::find(thisLayer.InfoNamesBegin(), thisLayer.InfoNamesEnd(), "NURBS_DEGREE_2");
                    if(it3 != thisLayer.InfoNamesEnd())
                    {
                        EntityInfoIterator it4 = (*it2)->find(*it3);
                        // check if this entity contain the info
                        if(it4 != (*it2)->end())
                            order_v = boost::get<int>(*it4);
                    }
                    it3 = std::find(thisLayer.InfoNamesBegin(), thisLayer.InfoNamesEnd(), "NURBS_DEGREE_3");
                    if(it3 != thisLayer.InfoNamesEnd())
                    {
                        EntityInfoIterator it4 = (*it2)->find(*it3);
                        // check if this entity contain the info
                        if(it4 != (*it2)->end())
                            order_w = boost::get<int>(*it4);
                    }

                    rOStream << order_u << " " << order_v << " " << order_w << std::endl;

                    //write the weights
                    EntityInfoIterator it4 = (*it2)->find("NURBS_WEIGHT");
                    if(it4 == (*it2)->end())
                        KRATOS_THROW_ERROR(std::runtime_error, "NURBS_WEIGHT is not provided for entity", (*it2)->Id())
                    rOStream << "        " << *it4 << std::endl;

                    //write the extraction operator
                    if((it4 = (*it2)->find("EXTRACTION_OPERATOR")) != (*it2)->end())
                    {
                        rOStream << "        Full" << std::endl;
                        rOStream << "        " << *it4 << std::endl;
                    }
                    else if((it4 = (*it2)->find("EXTRACTION_OPERATOR_MCSR")) != (*it2)->end())
                    {
                        rOStream << "        MCSR" << std::endl;
                        rOStream << "        " << *it4 << std::endl;
                    }
                    else if((it4 = (*it2)->find("EXTRACTION_OPERATOR_CSR")) != (*it2)->end())
                    {
//                        rOStream << "        CSR" << std::endl;
//                        rOStream << "        " << *it4 << std::endl;
                        KRATOS_THROW_ERROR(std::runtime_error, "EXTRACTION_OPERATOR_CSR is not implemented yet. Please wat.", "")
                    }
                    else
                        KRATOS_THROW_ERROR(std::runtime_error, "EXTRACTION_OPERATOR|EXTRACTION_OPERATOR_MCSR|EXTRACTION_OPERATOR_CSR is not provided for entity", (*it2)->Id())

                    rOStream << std::endl;
                }

                rOStream << "    ";
                End(rOStream, "IsogeometricBezierData");

                rOStream << "    ";
                Begin(rOStream, "ConditionsWithGeometry", layer_entity_name);
                for(Layer::EntitiesIteratorType it2 = thisLayer.EntitiesBegin(); it2 != thisLayer.EntitiesEnd(); ++it2)
                {
                    rOStream << "        " << (*it2)->Id() << " " << layer_prop_id << " " << (*it2)->Id();
                    for(std::size_t i = 0; i < (*it2)->size(); ++i)
                        rOStream << " " << (*it2)->operator[](i).Id();
                    rOStream << std::endl;
                }
                rOStream << "    ";
                End(rOStream, "ConditionsWithGeometry");

                End(rOStream, "BezierBlock");

                //write the entity info
                for(Layer::InfoNamesIteratorType it2 = thisLayer.InfoNamesBegin(); it2 != thisLayer.InfoNamesEnd(); ++it2)
                {
                    if(it2->find("EXTRACTION_OPERATOR") != std::string::npos
                        || it2->find("NURBS_WEIGHT") != std::string::npos
                        || it2->find("NURBS_DEGREE") != std::string::npos )
                    {
                        continue;
                    }

                    Begin(rOStream, "ConditionalData", *it2);
                    for(Layer::EntitiesIteratorType it3 = thisLayer.EntitiesBegin(); it3 != thisLayer.EntitiesEnd(); ++it3)
                    {
                        EntityInfoIterator it4 = (*it3)->find(*it2);
                        // check if this entity contain the info
                        if(it4 != (*it3)->end())
                        {
                            rOStream << (*it3)->Id() << " " << *it4 << std::endl;
                        }
                    }
                    End(rOStream, "ConditionalData");
                }

                std::cout << "Layer name " << thisLayer.Name() << " of type " << layer_entity_name << "(" << layer_entity_type << ") is written to Mdpa" << std::endl;
            }
            else if(layer_entity_type == std::string("null")) // null element/condition will not be accounted
            {
                continue;
            }
            else
            {
                std::stringstream ss;
                ss << "Unsupported entity type " << layer_entity_type << " of Layer " << thisLayer.Name();
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
        }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    LayerHandler& operator=(LayerHandler const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    LayerHandler(LayerHandler const& rOther)
    {
    }

    ///@}

}; // Class LayerHandler

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, LayerHandler& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const LayerHandler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_LAYER_HANDLER_H_INCLUDED
