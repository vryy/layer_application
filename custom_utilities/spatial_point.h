//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Sep 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_SPATIAL_POINT_H_INCLUDED )
#define KRATOS_LAYER_APP_SPATIAL_POINT_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "utilities/indexed_object.h"

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

template<class TDataType>
class SpatialPoint
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SpatialPoint);
    typedef std::size_t IndexType;

    // default constructor
    SpatialPoint()
    : mId(0), mLocalId(0), mLayerId(0), mX(0.0), mY(0.0), mZ(0.0)
    {}

    // constructor for pointer_vector_set which uses indexed_object
    SpatialPoint(IndexType NewId)
    : mId(NewId), mLocalId(0), mLayerId(0), mX(0.0), mY(0.0), mZ(0.0)
    {}

    // main constructor
    SpatialPoint(IndexType NewId, IndexType LocalId, IndexType LayerId, TDataType X, TDataType Y, TDataType Z)
    : mId(NewId), mLocalId(LocalId), mLayerId(LayerId), mX(X), mY(Y), mZ(Z)
    {}

    ~SpatialPoint() {}

    TDataType X() const {return mX;}
    TDataType Y() const {return mY;}
    TDataType Z() const {return mZ;}

    void SetId(IndexType Id) {mId = Id;}
    IndexType Id() const {return mId;}

    IndexType LocalId() const {return mLocalId;}
    IndexType LayerId() const {return mLayerId;}

    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SpatialPoint #" << Id();
    }

    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "(" << X() << ", " << Y() << ", " << Z() << ")";
    }

private:
    TDataType mX, mY, mZ;
    IndexType mId; // the global id of this node
    IndexType mLayerId; // the id of the layer containing this node
    IndexType mLocalId; // the unique local id of this node in the layer

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("Id", mId);
        rSerializer.save("LocalId", mLocalId);
        rSerializer.save("LayerId", mLayerId);
        rSerializer.save("X", mX);
        rSerializer.save("Y", mY);
        rSerializer.save("Z", mZ);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("Id", mId);
        rSerializer.load("LocalId", mLocalId);
        rSerializer.load("LayerId", mLayerId);
        rSerializer.load("X", mX);
        rSerializer.load("Y", mY);
        rSerializer.load("Z", mZ);
    }

}; // class SpatialPoint

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
inline std::ostream& operator <<(std::ostream& rOStream, const SpatialPoint<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_SPATIAL_POINT_H_INCLUDED defined
