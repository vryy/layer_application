//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_SPATIAL_KEY_H_INCLUDED )
#define KRATOS_LAYER_APP_SPATIAL_KEY_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"


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
 * This class defines a key in space. Used for spatial binning algorithm.
 */
class SpatialKey
{
    public:
        KRATOS_CLASS_POINTER_DEFINITION(SpatialKey);
        SpatialKey(int ix, int iy, int iz) : x(ix), y(iy), z(iz) {}
        ~SpatialKey() {}
        bool operator<(const SpatialKey& rOther) const
        {
            if(x == rOther.x)
            {
                if(y == rOther.y)
                {
                    return z < rOther.z;
                }
                else
                    return y < rOther.y;
            }
            else
                return x < rOther.x;
        }
        int kx() const {return x;}
        int ky() const {return y;}
        int kz() const {return z;}
    private:
        int x, y, z;
};

/// input stream function
// inline std::istream& operator >>(std::istream& rIStream, SpatialKey& rThis)
// {
//    return rIStream;
// }

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const SpatialKey& rThis)
{
    rOStream << '(' << rThis.kx() << ',' << rThis.ky() << ',' << rThis.kz() << ')';
    return rOStream;
}

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_SPATIAL_KEY_H_INCLUDED defined
