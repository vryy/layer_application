//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 31 Oct 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_LAYER_H_INCLUDED )
#define  KRATOS_LAYER_APP_LAYER_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <set>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/spatial_point.h"
#include "custom_utilities/entity.h"
#include "utilities/indexed_object.h"
#include "layer_application.h"

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

template<class TObjectType>
class LocalIdentFunction
{
public:
    typedef typename TObjectType::IndexType IndexType;
    typedef IndexType result_type;

    IndexType operator()(TObjectType const& rThisObject) const
    {
        return rThisObject.LocalId();
    }
};


/// Short class definition.
/*** Detail class definition.
    Layer implementation to hold nodes identification and entities's connectivities
 */
class Layer : public ParameterList<std::string>
{
public:

    ///@name Type Definitions
    ///@{

    typedef Entity EntityType;
    typedef EntityType::PointType PointType;
    typedef EntityType::IndexType IndexType;
    typedef ParameterList<std::string> BaseType;

//    typedef std::vector<PointType::Pointer> NodesContainerType;
//    typedef NodesContainerType::iterator NodesIteratorType;
    typedef PointerVectorSet<PointType, LocalIdentFunction<PointType> > NodesContainerType;
    typedef NodesContainerType::ptr_iterator NodesIteratorType;
    typedef NodesContainerType::ptr_const_iterator NodesConstIteratorType;

    typedef std::vector<EntityType::Pointer> EntitiesContainerType;
    typedef EntitiesContainerType::iterator EntitiesIteratorType;
    typedef std::set<std::string>::iterator InfoNamesIteratorType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Layer);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Layer(std::size_t Id, std::string Name) : mName(Name), mId(Id)
    {
        (*this)["LAYER_ENTITY_TYPE"] = std::string("element");
        (*this)["LAYER_ENTITY_NAME"] = std::string("Entity");
    }

    /// Destructor.
    virtual ~Layer()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    NodesContainerType& Nodes()
    {
        return mpNodes;
    }

    NodesIteratorType NodesBegin()
    {
        return mpNodes.ptr_begin();
    }

    NodesConstIteratorType NodesBegin() const
    {
        return mpNodes.ptr_begin();
    }

    NodesIteratorType NodesEnd()
    {
        return mpNodes.ptr_end();
    }

    NodesConstIteratorType NodesEnd() const
    {
        return mpNodes.ptr_end();
    }

    EntitiesIteratorType EntitiesBegin()
    {
        return mpEntities.begin();
    }

    EntitiesIteratorType EntitiesEnd()
    {
        return mpEntities.end();
    }

    EntitiesIteratorType find(IndexType id)
    {
        EntitiesIteratorType it = EntitiesBegin();
        while(it != EntitiesEnd())
        {
            if((*it)->Id() == id)
                break;
            else
                ++it;
        }
        return it;
    }

    InfoNamesIteratorType InfoNamesBegin()
    {
        return mInfoNames.begin();
    }

    InfoNamesIteratorType InfoNamesEnd()
    {
        return mInfoNames.end();
    }

    std::vector<std::string> GetTableNames() const
    {
        std::vector<std::string> table_names;
        for(std::map<std::string, std::vector<PointType::Pointer> >::const_iterator it = mpTables.begin(); it != mpTables.end(); ++it)
            table_names.push_back(it->first);
        return table_names;
    }

    std::vector<PointType::Pointer>& Table(std::string name)
    {
        std::map<std::string, std::vector<PointType::Pointer> >::iterator it = mpTables.find(name);
        if(it == mpTables.end())
        {
            std::stringstream ss;
            ss << "Table " << name << " does not exist in the layer";
            KRATOS_THROW_ERROR(std::runtime_error, ss.str(), "")
        }
        return it->second;
    }

    ///@}
    ///@name Operations
    ///@{

    void ClearNodes()
    {
        mpNodes.clear();
    }

    void ClearEntities()
    {
        mpEntities.clear();
    }

    void AddNode(PointType::Pointer pPoint)
    {
        mpNodes.push_back(pPoint);
    }

    void AddEntity(EntityType::Pointer pEntity)
    {
        mpEntities.push_back(pEntity);
    }

    void AddGroup(std::string Name)
    {
        mGroups.insert(Name);
    }

    void AddInfo(std::string Name)
    {
        mInfoNames.insert(Name);
    }

    void AddTable(std::string table_name, boost::python::list& pyListNodes)
    {
        std::vector<PointType::Pointer>& table_nodes = mpTables[table_name];
        typedef boost::python::stl_input_iterator<int> iterator_type;
        BOOST_FOREACH(const iterator_type::value_type& id,
                      std::make_pair(iterator_type(pyListNodes), // begin
                        iterator_type() ) ) // end
        {
            table_nodes.push_back(mpNodes(id));
        }
        std::cout << "Table " << table_name << " is added to layer " << mName << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    std::size_t NumberOfNodes() const
    {
        return mpNodes.size();
    }

    std::size_t NumberOfEntities() const
    {
        return mpEntities.size();
    }

    std::string Name() const
    {
        return mName;
    }

    std::size_t Id() const
    {
        return mId;
    }

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
        buffer << "Layer " << "'" << mName << "'";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "-->Nodes:";
        for(NodesConstIteratorType it = NodesBegin(); it != NodesEnd(); ++it)
            rOStream << " " << *(*it);
        rOStream << std::endl;
        rOStream << "-->Entities:" << std::endl;
        for(EntitiesContainerType::const_iterator it = mpEntities.begin(); it != mpEntities.end(); ++it)
        {
            rOStream << "   " << *(*it) << std::endl;
        }
        rOStream << "-->Entity Info:" << std::endl;
        for(EntitiesContainerType::const_iterator it = mpEntities.begin(); it != mpEntities.end(); ++it)
        {
            rOStream << "   " << (*it)->Id() << ":";
            rOStream << "   " << *(*it);
            rOStream << std::endl;
        }
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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

    std::string mName;
    std::size_t mId;
    NodesContainerType mpNodes;
    EntitiesContainerType mpEntities;
    std::set<std::string> mGroups; // list of group contain this layer
    std::set<std::string> mInfoNames; // list of info names associated with this layer
    std::map<std::string, std::vector<PointType::Pointer> > mpTables; // container of the subset of nodes in the layer (i.e to contain the boundary nodes)

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
    Layer& operator=(Layer const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    Layer(Layer const& rOther)
    {
    }

    ///@}

}; // Class Layer

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, Layer& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const Layer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_LAYER_H_INCLUDED
