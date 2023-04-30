//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 5 May 2016 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_ENTITY_H_INCLUDED )
#define KRATOS_LAYER_APP_ENTITY_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#ifdef SD_APP_FORWARD_COMPATIBILITY
#include "includes/indexed_object.h"
#else
#include "utilities/indexed_object.h"
#endif

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

class Entity : public ParameterList<std::string>, public IndexedObject
{
public:
    typedef ParameterList<std::string> BaseType;
    typedef SpatialPoint<double> PointType;
    typedef IndexedObject::IndexType IndexType;
    KRATOS_CLASS_POINTER_DEFINITION(Entity);

    Entity(IndexType NewId) : IndexedObject(NewId)
    {}

    ~Entity()
    {}

    void AddNode(PointType::Pointer pPoint)
    {
        mpConnectivities.push_back(pPoint);
    }

    std::size_t size() const
    {
        return mpConnectivities.size();
    }

    PointType& operator[] (std::size_t i)
    {
        return *mpConnectivities[i];
    }

    BaseType::DataType& operator[] (const std::string& str)
    {
        return BaseType::operator[](str);
    }

    bool operator< (const Entity& rOther) const
    {
        return Id() < rOther.Id();
    }

    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Entity #" << Id();
    }

    void PrintData(std::ostream& rOStream) const
    {
        for(std::size_t i = 0; i < this->size(); ++i)
            rOStream << " " << mpConnectivities[i]->Id();
    }

private:
    std::vector<PointType::Pointer> mpConnectivities;
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >>(std::istream& rIStream, SpatialPoint& rThis)
// {
//    return rIStream;
// }
//
/// output stream function
template<class TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const Entity& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_ENTITY_H_INCLUDED defined
