//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jul 2019 $
//   Revision:            $Revision: 1.0 $
//
//
// Project includes
#include "custom_utilities/model_part_utilities.h"


namespace Kratos
{

const int ModelPartUtilities::msT3Edges[][2] = { {0, 1}, {1, 2}, {2, 3} };
const int ModelPartUtilities::msQ4Edges[][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
const int ModelPartUtilities::msT4Edges[][2] = { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };
const int ModelPartUtilities::msH8Edges[][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7} };

}// namespace Kratos.


