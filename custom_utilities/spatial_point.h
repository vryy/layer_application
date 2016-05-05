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
class SpatialPoint : public IndexedObject
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SpatialPoint);
    typedef IndexedObject BaseType;
    typedef BaseType::IndexType IndexType;
    SpatialPoint() : BaseType(0), mX(0.0), mY(0.0), mZ(0.0) {}
    SpatialPoint(IndexType NewId) : BaseType(NewId), mX(0.0), mY(0.0), mZ(0.0) {}
    SpatialPoint(IndexType NewId, TDataType X, TDataType Y, TDataType Z) : BaseType(NewId), mX(X), mY(Y), mZ(Z) {}
    ~SpatialPoint() {}
    TDataType X() const {return mX;}
    TDataType Y() const {return mY;}
    TDataType Z() const {return mZ;}

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
