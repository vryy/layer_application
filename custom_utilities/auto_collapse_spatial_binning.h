//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 18 Sep 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED )
#define  KRATOS_LAYER_APP_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "spatial_point.h"
#include "spatial_key.h"


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
 * This utility class supports for spatial binning with auto collapsing functionality. This class used SpatialPoint as point data and will collapse the conincident point.
 */
class AutoCollapseSpatialBinning
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef SpatialPoint<double> SpatialPointType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(AutoCollapseSpatialBinning);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AutoCollapseSpatialBinning(double X0, double Y0, double Z0, double Dx, double Dy, double Dz, double tol)
    : mX0(X0), mY0(Y0), mZ0(Z0), mDx(Dx), mDy(Dy), mDz(Dz), mTol(tol)
    {
        mLastNode = 0;
    }

    /// Destructor.
    virtual ~AutoCollapseSpatialBinning()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function add node to the spatial binning and return the id of this node in the bin
     */
    IndexType AddNode(double X, double Y, double Z)
    {
        // find the cell containing point
        int ix = (int) floor((X - mX0) / mDx);
        int iy = (int) floor((Y - mY0) / mDy);
        int iz = (int) floor((Z - mZ0) / mDz);

        // check if the spatial key exist
        SpatialKey key(ix, iy, iz);
    //    std::cout << "ix = " << ix << ", iy = " << iy << ", iz = " << iz << std::endl;
        std::map<SpatialKey, std::vector<int> >::iterator it = mBin.find(key);
        if(it != mBin.end())
        {
            // check if node already exist in the cell
    //        std::cout << "check if node already exist in the cell" << std::endl;
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                double xi = mPointList[it->second[i] - 1]->X();
                double yi = mPointList[it->second[i] - 1]->Y();
                double zi = mPointList[it->second[i] - 1]->Z();
                double d = sqrt(pow(X - xi, 2) + pow(Y - yi, 2) + pow(Z - zi, 2));
                if(d < mTol)
                {
    //                std::cout << "node exist in cell, return its id" << std::endl;
                    return mPointList[it->second[i] - 1]->Id();
                }
            }
            // node does not exist in cell, insert the node into spatial bin
            mPointList.push_back(boost::shared_ptr<SpatialPointType>(new SpatialPointType(++mLastNode, 0, 0, X, Y, Z)));
            mBin[key].push_back(mLastNode);
    //        std::cout << "node does not exist in cell, insert the node into spatial bin, mLastNode = " << mLastNode << ", mBin[key].size() = " << mBin[key].size() << std::endl;
            return mLastNode;
        }
        else
        {
            // insert the node into spatial bin
            mPointList.push_back(SpatialPointType::Pointer(new SpatialPointType(++mLastNode, 0, 0, X, Y, Z)));
            mBin[key].push_back(mLastNode);
    //        std::cout << "insert the node into spatial bin, mLastNode = " << mLastNode << ", mBin[key].size() = " << mBin[key].size() << std::endl;
            return mLastNode;
        }
    }

    std::size_t NumberOfNodes() const {return mPointList.size();}
    double GetX(std::size_t id) const {return mPointList[id - 1]->X();}
    double GetY(std::size_t id) const {return mPointList[id - 1]->Y();}
    double GetZ(std::size_t id) const {return mPointList[id - 1]->Z();}

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
        buffer << "AutoCollapseSpatialBinning";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AutoCollapseSpatialBinning";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

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
    double mX0, mY0, mZ0, mDx, mDy, mDz;
    IndexType mLastNode;
    double mTol;

    std::map<SpatialKey, std::vector<int> > mBin;
    std::vector<boost::shared_ptr<SpatialPointType> > mPointList;

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
    AutoCollapseSpatialBinning& operator=(AutoCollapseSpatialBinning const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    AutoCollapseSpatialBinning(AutoCollapseSpatialBinning const& rOther)
    {
    }

    ///@}

}; // Class AutoCollapseSpatialBinning

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
//inline std::istream& operator >>(std::istream& rIStream, AutoCollapseSpatialBinning& rThis)
//{
//    return rIStream;
//}

///// output stream function
//inline std::ostream& operator <<(std::ostream& rOStream,
//        const AutoCollapseSpatialBinning& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED defined
