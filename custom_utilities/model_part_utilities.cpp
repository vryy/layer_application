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

void ModelPartUtilities::ExportNodalCoordinates(std::ostream& rOStream, ModelPart& r_model_part,
    const bool with_id, const bool with_x, const bool with_y, const bool with_z,
    const std::string& separator, const int precision)
{
    rOStream << std::setprecision(precision);

    for(ModelPart::NodeIterator i_node = r_model_part.NodesBegin();
            i_node != r_model_part.NodesEnd(); ++i_node)
    {
        bool init = false;

        if (with_id)
        {
            rOStream << i_node->Id();
            init = true;
        }

        if (with_x)
        {
            if (init)
                rOStream << separator;
            rOStream << i_node->X0();
            init = true;
        }

        if (with_y)
        {
            if (init)
                rOStream << separator;
            rOStream << i_node->Y0();
        }

        if (with_z)
        {
            if (init)
                rOStream << separator;
            rOStream << i_node->Z0();
        }

        rOStream << std::endl;
    }
}

void ModelPartUtilities::ExportNodalCoordinatesToGiD(std::ostream& rOStream, ModelPart& r_model_part,
    const int precision)
{
    rOStream << std::setprecision(precision);

    for(ModelPart::NodeIterator i_node = r_model_part.NodesBegin();
            i_node != r_model_part.NodesEnd(); ++i_node)
    {
        rOStream << "Mescape Geometry Create Point "
                 << i_node->X0() << " " << i_node->Y0() << " " << i_node->Z0()
                 << std::endl;
    }
}

void ModelPartUtilities::ExportEdgeInformation(std::ostream& rOStream, ModelPart::ElementsContainerType& rpElements,
    const std::string& separator)
{
    typedef std::pair<std::size_t, std::size_t> edge_t;
    typedef std::set<edge_t> edge_container_t;

    edge_container_t edges;
    ExtractEdgeInformation(edges, rpElements);

    for (edge_container_t::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        rOStream << it->first << separator << it->second << std::endl;
    }
}

void ModelPartUtilities::ExportEdgeInformationToGiD(std::ostream& rOStream, ModelPart::ElementsContainerType& rpElements)
{
    typedef std::pair<std::size_t, std::size_t> edge_t;
    typedef std::set<edge_t> edge_container_t;

    edge_container_t edges;
    ExtractEdgeInformation(edges, rpElements);

    for (edge_container_t::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        rOStream << "Mescape Geometry Create Line Join " << it->first << " " << it->second << std::endl;
    }
}

}// namespace Kratos.
