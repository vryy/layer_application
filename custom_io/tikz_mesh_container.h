//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Hoang-Giang Bui $
//   Date:                $Date: 8 Sep 2016 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_TIKZ_MESH_CONTAINER_H_INCLUDED)
#define  KRATOS_TIKZ_MESH_CONTAINER_H_INCLUDED
// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"
#include "includes/deprecated_variables.h"


namespace Kratos
{

/**
 * Type definitions
 */
typedef ModelPart::ElementsContainerType ElementsArrayType;
typedef ModelPart::NodesContainerType NodesArrayType;
typedef ModelPart::ConditionsContainerType ConditionsArrayType;
typedef GeometryData::IntegrationMethod IntegrationMethodType;
typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

/**
 * Auxiliary class to store meshes of different element types and to
 * write these meshes to an output file
 */
class TikzMeshContainer
{
public:

    class Line
    {
        public:
            Line(const std::size_t& N1, const std::size_t& N2)
            {
                if(N1 > N2)
                {
                    mN1 = N2;
                    mN2 = N1;
                }
                else
                {
                    mN1 = N1;
                    mN2 = N2;
                }
            }
            const std::size_t& N1() const {return mN1;}
            const std::size_t& N2() const {return mN2;}
            bool operator<(const Line& rOther) const
            {
                if(mN1 == rOther.mN1)
                    return mN2 < rOther.mN2;
                else
                    return mN1 < rOther.mN1;
            }
        private:
            std::size_t mN1, mN2;
    };

    KRATOS_CLASS_POINTER_DEFINITION(TikzMeshContainer);

    ///Constructor
    TikzMeshContainer ( GeometryData::KratosGeometryType geometryType, const char* mesh_title )
    : mMeshTitle (mesh_title)
    {
        mGeometryType = geometryType;
    }

    bool AddElement ( const ModelPart::ElementsContainerType::iterator pElemIt )
    {
        KRATOS_TRY

        if ( pElemIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
//            std::cout << "element " << pElemIt->Id() << " of geometryType "
//                      << mGeometryType << " is added to " << mMeshTitle << std::endl;
            mMeshElements.push_back ( * (pElemIt.base() ) );
            Geometry<Node<3> >&geom = pElemIt->GetGeometry();
            for ( Element::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;

        KRATOS_CATCH ("")
    }

    bool AddCondition (const ModelPart::ConditionsContainerType::iterator pCondIt)
    {
        KRATOS_TRY

        if ( pCondIt->GetGeometry().GetGeometryType() == mGeometryType )
        {
//            std::cout << "condition " << pCondIt->Id() << " of geometryType "
//                      << mGeometryType << " is added to " << mMeshTitle << std::endl;
            mMeshConditions.push_back ( * (pCondIt.base() ) );
            Geometry<Node<3> >&geom = pCondIt->GetGeometry();
            for ( Condition::GeometryType::iterator it = geom.begin(); it != geom.end(); it++)
            {
                mMeshNodes.push_back ( * (it.base() ) );
            }
            return true;
        }
        else
            return false;

        KRATOS_CATCH ("")
    }

    void FinalizeMeshCreation()
    {
        if ( mMeshElements.size() != 0 )
        {
            mMeshNodes.Unique();
        }
        if ( mMeshConditions.size() != 0 )
        {
            mMeshNodes.Unique();
        }
    }

    void WriteMesh(std::ostream& rMeshFile,
            bool deformed, bool write_id,
            std::string node_style, double node_size,
            std::string element_style, std::string condition_style)
    {
        KRATOS_TRY

        bool nodes_written = false;
        if ( mMeshElements.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> elements_per_layer (max_id+1,0);

            //fill layer list
            for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                    it != mMeshElements.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                elements_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing elements" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < elements_per_layer.size(); current_layer++)
            {
                if (elements_per_layer[current_layer] > 0)
                {
                    //printing nodes
                    if(nodes_written == false)
                    {
                        for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
                        {
                            if ( deformed )
                                rMeshFile << "\\coordinate (P" << it->Id() << ") at (" << it->X() << ", " << it->Y() << ", " << it->Z() << ");" << std::endl;
                            else
                                rMeshFile << "\\coordinate (P" << it->Id() << ") at (" << it->X0() << ", " << it->Y0() << ", " << it->Z0() << ");" << std::endl;
                            if(write_id)
                                rMeshFile << "\\draw[" << node_style << "] (P" << it->Id() << ") circle (" << node_size << "em) node[above right] {" << it->Id() << "};" << std::endl;
                            else
                                rMeshFile << "\\draw[" << node_style << "] (P" << it->Id() << ") circle (" << node_size << "em) node[above right] {};" << std::endl;
                        }

                        nodes_written = true;
                    }

                    //printing elements
                    for ( ModelPart::ElementsContainerType::iterator it = mMeshElements.begin();
                            it != mMeshElements.end(); ++it )
                    {
                        if ( it->Has ( IS_INACTIVE ) )
                        {
                            if ( ! it->GetValue ( IS_INACTIVE )  && (it)->GetProperties().Id()==current_layer )
                            {
                                DrawGeometry(rMeshFile, (it)->GetGeometry(), element_style);
                                if(write_id)
                                {
                                    Element::GeometryType::PointType C = it->GetGeometry().Center();
                                    rMeshFile << "\\draw (" << C.X() << "," << C.Y() << "," << C.Z() << ") node[above] {" << it->Id() << "};" << std::endl;
                                }
                            }
                        }
                        else
                        {
                            if ((it)->GetProperties().Id()==current_layer)
                            {
                                DrawGeometry(rMeshFile, (it)->GetGeometry(), element_style);
                                if(write_id)
                                {
                                    Element::GeometryType::PointType C = it->GetGeometry().Center();
                                    rMeshFile << "\\draw (" << C.X() << "," << C.Y() << "," << C.Z() << ") node[above] {" << it->Id() << "};" << std::endl;
                                }
                            }
                        }
                    }
                }
            }
            //std::cout << "end printing elements" <<std::endl;
        }

        if ( mMeshConditions.size() != 0 )
        {
            //compute number of layers
            int max_id = 0;
            for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                if (max_id < prop_id) max_id = prop_id;
            }
            if (max_id > 10000)
                std::cout<< "a property Id > 10000 found. Are u sure you need so many properties?" << std::endl;
            std::vector<int> conditions_per_layer (max_id+1,0);
            //fill layer list
            for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin();
                    it != mMeshConditions.end(); ++it )
            {
                int prop_id = (it)->GetProperties().Id();
                conditions_per_layer[prop_id] += 1;
            }
            //std::cout << "start printing conditions" <<std::endl;
            for (unsigned int current_layer = 0; current_layer < conditions_per_layer.size(); current_layer++)
            {
                if (conditions_per_layer[current_layer] > 0)
                {
                    std::stringstream current_layer_name (std::stringstream::in | std::stringstream::out);
                    current_layer_name << mMeshTitle << "_" << current_layer ;

                    //printing nodes
                    if(nodes_written == false)
                    {
                        for ( ModelPart::NodesContainerType::iterator it = mMeshNodes.begin();
                                it != mMeshNodes.end(); ++it )
                        {
                            if ( deformed )
                                rMeshFile << "\\coordinate (P" << it->Id() << ") at (" << it->X() << ", " << it->Y() << ", " << it->Z() << ");" << std::endl;
                            else
                                rMeshFile << "\\coordinate (P" << it->Id() << ") at (" << it->X0() << ", " << it->Y0() << ", " << it->Z0() << ");" << std::endl;
                            if(write_id)
                                rMeshFile << "\\draw[" << node_style << "] (P" << it->Id() << ") circle (" << node_size << "em) node[above right] {" << it->Id() << "};" << std::endl;
                            else
                                rMeshFile << "\\draw[" << node_style << "] (P" << it->Id() << ") circle (" << node_size << "em) node[above right] {};" << std::endl;
                        }

                        nodes_written = true;
                    }

                    //printing conditions
                    for ( ModelPart::ConditionsContainerType::iterator it = mMeshConditions.begin(  );
                            it != mMeshConditions.end(); ++it )
                    {
                        if ( it->Has ( IS_INACTIVE ) )
                        {
                            if ( ! it->GetValue ( IS_INACTIVE )  && (it)->GetProperties().Id()==current_layer )
                            {
                                DrawGeometry(rMeshFile, (it)->GetGeometry(), condition_style);
                                if(write_id)
                                {
                                    Element::GeometryType::PointType C = (it)->GetGeometry().Center();
                                    rMeshFile << "\\draw (" << C.X() << "," << C.Y() << "," << C.Z() << ") node[above] {" << it->Id() << "};" << std::endl;
                                }
                            }
                        }
                        else
                        {
                            if ((it)->GetProperties().Id()==current_layer )
                            {
                                DrawGeometry(rMeshFile, (it)->GetGeometry(), condition_style);
                                if(write_id)
                                {
                                    Element::GeometryType::PointType C = (it)->GetGeometry().Center();
                                    rMeshFile << "\\draw (" << C.X() << "," << C.Y() << "," << C.Z() << ") node[above] {" << it->Id() << "};" << std::endl;
                                }
                            }

                        }
                    }
                }
            }
            //std::cout << "end printing conditions" <<std::endl;
        }

        KRATOS_CATCH ("")
    }

    void Reset()
    {
        mMeshNodes.clear();
        mMeshElements.clear();
        mMeshConditions.clear();
    }

    ModelPart::NodesContainerType GetMeshNodes()
    {
        return mMeshNodes;
    }

protected:
    ///member variables
    GeometryData::KratosGeometryType mGeometryType;
    ModelPart::NodesContainerType mMeshNodes;
    ModelPart::ElementsContainerType mMeshElements;
    ModelPart::ConditionsContainerType mMeshConditions;
    const char* mMeshTitle;

private:
    void AddLines(std::set<Line>& rlines, std::size_t* nodes, GeometryData::KratosGeometryType GeometryType)
    {
        if(GeometryType == GeometryData::Kratos_Hexahedra3D8)
        {
            rlines.insert(Line(nodes[0], nodes[1]));
            rlines.insert(Line(nodes[1], nodes[2]));
            rlines.insert(Line(nodes[2], nodes[3]));
            rlines.insert(Line(nodes[3], nodes[0]));
            rlines.insert(Line(nodes[4], nodes[5]));
            rlines.insert(Line(nodes[5], nodes[6]));
            rlines.insert(Line(nodes[6], nodes[7]));
            rlines.insert(Line(nodes[7], nodes[4]));
            rlines.insert(Line(nodes[0], nodes[4]));
            rlines.insert(Line(nodes[1], nodes[5]));
            rlines.insert(Line(nodes[2], nodes[6]));
            rlines.insert(Line(nodes[3], nodes[7]));
        }
        else if(GeometryType == GeometryData::Kratos_Quadrilateral3D4)
        {
            rlines.insert(Line(nodes[0], nodes[1]));
            rlines.insert(Line(nodes[1], nodes[2]));
            rlines.insert(Line(nodes[2], nodes[3]));
            rlines.insert(Line(nodes[3], nodes[0]));
        }
    }

    void DrawGeometry(std::ostream& rOStream, Element::GeometryType& rGeometry, std::string style)
    {
        rOStream << "\\draw[" << style << "]";

        if(     rGeometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4
            ||  rGeometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4
            ||  rGeometry.GetGeometryType() == GeometryData::Kratos_Triangle2D3
            ||  rGeometry.GetGeometryType() == GeometryData::Kratos_Triangle3D3
        )
        {
            rOStream << " (P" << rGeometry[0].Id() << ")";
            for(unsigned int i = 1; i < rGeometry.size(); ++i)
                rOStream << " -- (P" << rGeometry[i].Id() << ")";
        }
        else if(rGeometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral3D8
            ||  rGeometry.GetGeometryType() == GeometryData::Kratos_Quadrilateral3D9)
        {
            rOStream << " (P" << rGeometry[0].Id() << ")";
            rOStream << " -- (P" << rGeometry[4].Id() << ")";
            rOStream << " -- (P" << rGeometry[1].Id() << ")";
            rOStream << " -- (P" << rGeometry[5].Id() << ")";
            rOStream << " -- (P" << rGeometry[2].Id() << ")";
            rOStream << " -- (P" << rGeometry[6].Id() << ")";
            rOStream << " -- (P" << rGeometry[3].Id() << ")";
            rOStream << " -- (P" << rGeometry[7].Id() << ")";
        }

        rOStream << " --cycle;" << std::endl;
    }

};//class TikzMeshContainer

}// namespace Kratos.

#endif // KRATOS_TIKZ_MESH_CONTAINER_H_INCLUDED defined

