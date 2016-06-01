//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 3 May 2016 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LAYER_MULTIPATCH_IGA_LAYER_HANDLER_H_INCLUDED )
#define  KRATOS_LAYER_APP_MULTIPATCH_LAYER_HANDLER_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <ctime>

// External includes
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "custom_utilities/parameter_list.h"
#include "custom_utilities/layer.h"
#include "custom_utilities/mdpa_writer.h"


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

/// Short class definition.
/*** Detail class definition.
    Layer handler keeps a global node information and handle multipatch structure properly
 */
class MultipatchLayerHandler : public LayerHandler
{
public:

    ///@name Type Definitions
    ///@{

    typedef LayerHandler BaseType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::PointsContainerType PointsContainerType;
    typedef typename BaseType::LayersContainerType LayersContainerType;

    typedef std::pair<std::string, std::string> pair_t;
    typedef std::pair<std::size_t, std::size_t> point_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(MultipatchLayerHandler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultipatchLayerHandler() : BaseType(), mLastBucket(0)
    {
    }

    /// Destructor.
    virtual ~MultipatchLayerHandler()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void AddLayerConnection(std::string layer1_name, std::string layer1_boundary,
            std::string layer2_name, std::string layer2_boundary)
    {
        const double TOL = 1.0e-10;

        std::vector<Layer::PointType::Pointer>& layer1_nodes = mpLayers[layer1_name]->Table(layer1_boundary);
        std::vector<Layer::PointType::Pointer>& layer2_nodes = mpLayers[layer2_name]->Table(layer2_boundary);

        /* check the consistency between two boundaries */
        // first count the number of nodes on the boundary
        std::size_t num_nodes_1 = layer1_nodes.size();
        std::size_t num_nodes_2 = layer2_nodes.size();
        if(num_nodes_1 != num_nodes_2)
        {
            KRATOS_WATCH(num_nodes_1)
            KRATOS_WATCH(num_nodes_2)
            KRATOS_THROW_ERROR(std::logic_error, "The number of nodes on the coupling boundary are not the same", "")
        }

        // check the coordinates of each node one by one
        for(std::size_t i = 0; i < num_nodes_1; ++i)
        {
            Layer::PointType::Pointer p1 = layer1_nodes[i];
            Layer::PointType::Pointer p2 = layer2_nodes[i];

            double diff = sqrt(pow(p1->X() - p2->X(), 2) + pow(p1->Y() - p2->Y(), 2) + pow(p1->Z() - p2->Z(), 2));
            if(diff > TOL)
            {
                KRATOS_WATCH(diff)
                KRATOS_WATCH(*p1)
                KRATOS_WATCH(*p2)
                KRATOS_THROW_ERROR(std::logic_error, "The node coordinates are different", "")
            }

            // now we make the indicator of the two node the same
            point_t pnt1(p1->LayerId(), p1->LocalId());
            point_t pnt2(p2->LayerId(), p2->LocalId());

            if(mPointMap.find(pnt1) == mPointMap.end())
            {
                if(mPointMap.find(pnt2) == mPointMap.end())
                {
                    mPointMap[pnt1] = mLastBucket;
                    mPointMap[pnt2] = mLastBucket;
                    ++mLastBucket;
                }
                else
                {
                    mPointMap[pnt1] = mPointMap[pnt2];
                }
            }
            else
            {
                if(mPointMap.find(pnt2) == mPointMap.end())
                {
                    mPointMap[pnt2] = mPointMap[pnt1];
                }
                else
                {
                    std::size_t OldBucket = mPointMap[pnt2];
                    std::size_t NewBucket = mPointMap[pnt1];

                    for(std::map<point_t, std::size_t>::iterator it = mPointMap.begin(); it != mPointMap.end(); ++it)
                    {
                        if(it->second == OldBucket)
                            it->second = NewBucket;
                    }
                }
            }
        }

        std::cout << "Layer connection (" << layer1_name << ", " << layer1_boundary
                  << ")-(" << layer2_name << ", " << layer2_boundary << ") is added" << std::endl;
    }


    void FinalizeLayerConnection()
    {
//        for(std::map<point_t, std::size_t>::iterator it = mPointMap.begin(); it != mPointMap.end(); ++it)
//        {
//            Layer::Pointer& layer = this->operator[](it->first.first);
//            Layer::PointType& p = layer->Nodes()[it->first.second];

//            std::cout << "point " << p.Id() << " bucket " << it->second << std::endl;
//        }

        // extract the pairs with the same flag
        std::map<std::size_t, std::set<point_t> > point_set;
        for(std::map<point_t, std::size_t>::iterator it = mPointMap.begin(); it != mPointMap.end(); ++it)
        {
            point_set[it->second].insert(it->first);
        }

//        for(std::map<std::size_t, std::set<point_t> >::iterator it = point_set.begin(); it != point_set.end(); ++it)
//        {
//            std::cout << "bucket " << it->first << ":";
//            for(std::set<point_t>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
//            {
//                Layer::Pointer& layer = this->operator[](it2->first);
//                Layer::PointType& p = layer->Nodes()[it2->second];
//                std::cout << " " << p.Id();
//            }
//            std::cout << std::endl;
//        }

        // for each set of pair, try to equalize the id of each nodes in pair
        for(std::map<std::size_t, std::set<point_t> >::iterator it = point_set.begin(); it != point_set.end(); ++it)
        {
            std::vector<point_t> points(it->second.begin(), it->second.end());

            if(points.size() > 1)
            {
                Layer::Pointer& layer1 = this->operator[](points[0].first);
                Layer::PointType& p1 = layer1->Nodes()[points[0].second];

                for(std::size_t i = 1; i < points.size(); ++i)
                {
                    Layer::Pointer& layeri = this->operator[](points[i].first);
                    Layer::PointType& pi = layeri->Nodes()[points[i].second];

                    pi.SetId(p1.Id());
                }
            }
        }
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
        buffer << "MultipatchLayerHandler";
        return buffer.str();
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

    std::size_t mLastBucket;
    std::map<point_t, std::size_t> mPointMap;

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
    MultipatchLayerHandler& operator=(MultipatchLayerHandler const& rOther)
    {
        return *this;
    }

    /// Copy constructor.
    MultipatchLayerHandler(MultipatchLayerHandler const& rOther)
    {
    }

    ///@}

}; // Class MultipatchLayerHandler

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >>(std::istream& rIStream, MultipatchLayerHandler& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream, const MultipatchLayerHandler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@} addtogroup block

}// namespace Kratos.

#endif // KRATOS_LAYER_APP_MULTIPATCH_LAYER_HANDLER_H_INCLUDED
