//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 8 Nov 2014 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_APP_SPATIAL_GRID_NODAL_BINNING_H_INCLUDED )
#define KRATOS_LAYER_APP_SPATIAL_GRID_NODAL_BINNING_H_INCLUDED

// System includes
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
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
 * This utility class supports for spatial binning into a structured grid. This class use Node as point data and will not collapse the conincident node.
 */
class SpatialGridNodalBinning
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t IndexType;
    typedef std::map<SpatialKey, std::set<IndexType> > BinType;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(SpatialGridNodalBinning);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SpatialGridNodalBinning(double X0, double Y0, double Z0, double Dx, double Dy, double Dz, double tol)
    : mX0(X0), mY0(Y0), mZ0(Z0), mDx(Dx), mDy(Dy), mDz(Dz), mTol(tol)
    {
    }

    /// Destructor.
    virtual ~SpatialGridNodalBinning()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function add node to the spatial binning
     */
    void AddNodes(ModelPart& r_model_part)
    {
        for(typename ModelPart::NodeIterator it = r_model_part.NodesBegin();
            it != r_model_part.NodesEnd(); ++it)
        {
            // find the cell containing point
            int ix = (int) floor((it->X() - mX0) / mDx);
            int iy = (int) floor((it->Y() - mY0) / mDy);
            int iz = (int) floor((it->Z() - mZ0) / mDz);

            // create the spatial key
            SpatialKey key(ix, iy, iz);

            // insert the node into spatial bin
            mBin[key].insert(it->Id());
        }
    }

    /**
     * Get the list of neighbour nodes of node Id within radius r
     */
    std::set<IndexType> GetNeighbourNodes(ModelPart& r_model_part, const IndexType& id, const double& r)
    {
        // find the cell containing point
        double x = r_model_part.Nodes()[id].X();
        double y = r_model_part.Nodes()[id].Y();
        double z = r_model_part.Nodes()[id].Z();
        int ix = (int) floor((x - mX0) / mDx);
        int iy = (int) floor((y - mY0) / mDy);
        int iz = (int) floor((z - mZ0) / mDz);

        // compute the cell span
        int d_ix = (int) ceil(r / mDx);
        int d_iy = (int) ceil(r / mDy);
        int d_iz = (int) ceil(r / mDz);

        // iterate through cell span to find neighbours
        std::set<IndexType> Neighbours;
        for(int i = ix - d_ix; i <= ix + d_ix; ++i)
        {
            for(int j = iy - d_iy; j <= iy + d_iy; ++j)
            {
                for(int k = iz - d_iz; k <= iz + d_iz; ++k)
                {
                    SpatialKey key(i, j, k);
                    BinType::iterator it = mBin.find(key);
                    if(it != mBin.end())
                    {
                        std::set<IndexType> Nodes = mBin[key];
                        for(std::set<IndexType>::iterator it2 = Nodes.begin(); it2 != Nodes.end(); ++it2)
                        {
                            double dx = r_model_part.Nodes()[*it2].X() - x;
                            double dy = r_model_part.Nodes()[*it2].Y() - y;
                            double dz = r_model_part.Nodes()[*it2].Z() - z;
                            double dist = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
                            if(dist < r + mTol)
                                Neighbours.insert(*it2);
                        }
                    }
                }
            }
        }

        return Neighbours;
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
        buffer << "SpatialGridNodalBinning";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SpatialGridNodalBinning";
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

    ///@}
    ///@name Member Variables
    ///@{

    double mX0, mY0, mZ0, mDx, mDy, mDz;
    double mTol;

    BinType mBin;

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
    SpatialGridNodalBinning& operator=(SpatialGridNodalBinning const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    SpatialGridNodalBinning(SpatialGridNodalBinning const& rOther)
    {
    }

    ///@}

}; // Class SpatialGridNodalBinning

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
//inline std::istream& operator >>(std::istream& rIStream, SpatialGridNodalBinning& rThis)
//{
//    return rIStream;
//}

///// output stream function
//inline std::ostream& operator <<(std::ostream& rOStream,
//        const SpatialGridNodalBinning& rThis)
//{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

//    return rOStream;
//}
///@}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_SPATIAL_GRID_NODAL_BINNING_H_INCLUDED defined

