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
template<typename TIndexType = std::size_t, typename TDataType = double>
class AutoCollapseSpatialBinning
{
public:
    ///@name Type Definitions
    ///@{

    typedef SpatialPoint<TDataType> SpatialPointType;
    typedef std::map<SpatialKey, std::vector<TIndexType> > BinType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(AutoCollapseSpatialBinning);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AutoCollapseSpatialBinning(TDataType X0, TDataType Y0, TDataType Z0, TDataType Dx, TDataType Dy, TDataType Dz, TDataType tol)
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

    /// Set the binning distance
    void SetDistance(double dx, double dy, double dz)
    {
        mDx = dx;
        mDy = dy;
        mDz = dz;
    }

    /// Set the binning tolerance
    void SetTolerance(double tol)
    {
        mTol = tol;
    }

    /**
     * Add node to the spatial binning and return the id of this node in the bin
     */
    TIndexType AddNode(const TDataType X, const TDataType Y, const TDataType Z)
    {
        // find the cell containing point
        int ix = (int) floor((X - mX0) / mDx);
        int iy = (int) floor((Y - mY0) / mDy);
        int iz = (int) floor((Z - mZ0) / mDz);

        // check if the spatial key exist
        SpatialKey key(ix, iy, iz);

        typename BinType::iterator it = mBin.find(key);
        if(it != mBin.end())
        {
            // check if node already exist in the cell
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                const auto p = mPointList[it->second[i] - 1];
                TDataType xi = p->X();
                TDataType yi = p->Y();
                TDataType zi = p->Z();
                TDataType d = sqrt(pow(X - xi, 2) + pow(Y - yi, 2) + pow(Z - zi, 2));
                if(d < mTol)
                {
                    // node exist in cell, return its id
                    return p->Id();
                }
            }
            // node does not exist in cell, insert the node into spatial bin
            mPointList.push_back(typename SpatialPointType::Pointer(new SpatialPointType(++mLastNode, 0, 0, X, Y, Z)));
            mBin[key].push_back(mLastNode);
            return mLastNode;
        }
        else
        {
            // insert the node into spatial bin
            mPointList.push_back(typename SpatialPointType::Pointer(new SpatialPointType(++mLastNode, 0, 0, X, Y, Z)));
            mBin[key].push_back(mLastNode);
            return mLastNode;
        }
    }

    /**
     * Check if a node exists in bin
     */
    bool HasNode(const TDataType X, const TDataType Y, const TDataType Z) const
    {
        TIndexType dummy;
        return HasNode(X, Y, Z, dummy);
    }

    /**
     * Check if a node exists in bin. If existed, then return its index
     */
    bool HasNode(const TDataType X, const TDataType Y, const TDataType Z, TIndexType& node_index) const
    {
        // find the cell containing point
        int ix = (int) floor((X - mX0) / mDx);
        int iy = (int) floor((Y - mY0) / mDy);
        int iz = (int) floor((Z - mZ0) / mDz);

        // check if the spatial key exist
        SpatialKey key(ix, iy, iz);

        typename BinType::const_iterator it = mBin.find(key);
        if(it != mBin.end())
        {
            // check if node already exist in the cell
            for(std::size_t i = 0; i < it->second.size(); ++i)
            {
                const auto p = mPointList[it->second[i] - 1];
                TDataType xi = p->X();
                TDataType yi = p->Y();
                TDataType zi = p->Z();
                TDataType d = sqrt(pow(X - xi, 2) + pow(Y - yi, 2) + pow(Z - zi, 2));
                if(d < mTol)
                {
                    node_index = p->Id();
                    return true;
                }
            }
        }

        return false;
    }

    /**
     * Return the id of this node in the bin and throw error if that node does not exist
     */
    TIndexType GetNode(const TDataType X, const TDataType Y, const TDataType Z) const
    {
        TIndexType node_index;

        if (this->HasNode(X, Y, Z, node_index))
        {
            return node_index;
        }

        KRATOS_ERROR << "Node " << X << "," << Y << "," << Z << " does not exist";
        return -1; // silence the compiler
    }

    std::size_t NumberOfNodes() const {return mPointList.size();}
    TDataType GetX(std::size_t id) const {return mPointList[id - 1]->X();}
    TDataType GetY(std::size_t id) const {return mPointList[id - 1]->Y();}
    TDataType GetZ(std::size_t id) const {return mPointList[id - 1]->Z();}

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
    TDataType mX0, mY0, mZ0, mDx, mDy, mDz;
    TIndexType mLastNode;
    TDataType mTol;

    BinType mBin;
    std::vector<typename SpatialPointType::Pointer> mPointList;

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
template<typename TIndexType, typename TDataType>
inline std::istream& operator >>(std::istream& rIStream, AutoCollapseSpatialBinning<TIndexType, TDataType>& rThis)
{
   return rIStream;
}

/// output stream function
template<typename TIndexType, typename TDataType>
inline std::ostream& operator <<(std::ostream& rOStream, const AutoCollapseSpatialBinning<TIndexType, TDataType>& rThis)
{
   rThis.PrintInfo(rOStream);
   rOStream << std::endl;
   rThis.PrintData(rOStream);

   return rOStream;
}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_AUTO_COLLAPSE_SPATIAL_BINNING_H_INCLUDED defined
