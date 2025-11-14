//
//   Project Name:        KratosLayerApplication
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 14 Nov 2025 $
//
//

#if !defined(KRATOS_LAYER_APP_POINT_CLOUD_UTILITY_H_INCLUDED )
#define  KRATOS_LAYER_APP_POINT_CLOUD_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"


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

/**
 *  Utilities for operations on point cloud
 */
class PointCloudUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointCloudUtility);

    typedef KRATOS_INDEX_TYPE IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointCloudUtility()
    {
    }

    /// Destructor.
    virtual ~PointCloudUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Read from xyz file to point cloud
    /// Note:   - the line starting with # will be skipped
    ///         - the coordinate line should be id x y z
    template<typename TVectorType>
    static std::vector<TVectorType> ReadXYZ(const std::string& fileName)
    {
        std::vector<TVectorType> points;

        std::ifstream file(fileName);
        if (!file.is_open())
        {
            KRATOS_ERROR << "Could not open file " << fileName;
            return points;
        }

        std::string line;
        IndexType i;
        typename TVectorType::value_type x, y, z;

        while (std::getline(file, line))
        {
            // Trim leading spaces (optional)
            if (line.empty()) continue;
            if (line[0] == '#') continue;  // Skip comment lines

            std::istringstream iss(line);
            if (iss >> i >> x >> y >> z)
            {
                TVectorType p;
                p[0] = x;
                p[1] = y;
                p[2] = z;
                points.push_back(p);
            }
        }

        file.close();
        return points;
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
        buffer << "PointCloudUtility";
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
    PointCloudUtility& operator=(PointCloudUtility const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    PointCloudUtility(PointCloudUtility const& rOther)
    {
    }

    ///@}

}; // Class PointCloudUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, PointCloudUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const PointCloudUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_LAYER_APP_POINT_CLOUD_UTILITY_H_INCLUDED
