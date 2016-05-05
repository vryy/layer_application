//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 31 Oct 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_GROUP_H_INCLUDED )
#define  KRATOS_LAYER_APP_GROUP_H_INCLUDED

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
#include "custom_utilities/layer.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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
    Group implementation to hold layers and to handle collapse properly
 */
class Group
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::set<std::string> LayerContainerType;
    typedef std::set<std::string>::iterator LayerContainerIteratorType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(Group);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Group(std::string Name)
    {
        mName = Name;
    }

    /// Destructor.
    virtual ~Group()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    void AddLayer(Layer::Pointer& pLayer)
    {
        mLayers.insert(pLayer->Name());
        pLayer->AddGroup(mName);
    }
    
    ///@}
    ///@name Access
    ///@{

    LayerContainerIteratorType LayerBegin()
    {
        return mLayers.begin();
    }

    LayerContainerIteratorType LayerEnd()
    {
        return mLayers.end();
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
        buffer << "Group " << "'" << mName << "'";
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
        std::cout << "-->Group Layer:";
        for(std::set<std::string>::const_iterator it = mLayers.begin(); it != mLayers.end(); ++it)
            std::cout << " " << "'" << (*it) << "'";
        std::cout << std::endl;
    }

    void Print()
    {
        PrintInfo(std::cout);
        std::cout << std::endl;
        PrintData(std::cout);
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
    LayerContainerType mLayers;

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
    Group& operator=(Group const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    Group(Group const& rOther)
    {
    }

    ///@}

}; // Class Group

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, Group& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const Group& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_GROUP_H_INCLUDED

