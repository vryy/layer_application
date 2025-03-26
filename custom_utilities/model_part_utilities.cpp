//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jul 2019 $
//   Revision:            $Revision: 1.0 $
//
//
// Project includes
#include "custom_utilities/gidpost_binary_reader.h"
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

void ModelPartUtilities::GiDPostBin2ModelPart(const std::string& fileName, ModelPart& r_model_part,
        const Parameters& mesh_info)
{
    GiDPostBinaryReader reader(fileName);
    GiDPost2ModelPart(reader, r_model_part, mesh_info);
}

void ModelPartUtilities::GiDPost2ModelPart(GiDPostReader& reader, ModelPart& r_model_part,
        const Parameters& mesh_info)
{
    typedef ModelPart::IndexType IndexType;

    std::vector<std::string> mesh_names = reader.GetMeshesName();

    int echo_level = mesh_info.Has("echo_level") ? mesh_info.GetValue("echo_level").GetInt() : 0;

    for (const auto& name : mesh_names)
    {
        std::map<int, std::vector<double> > coordinates;
        std::map<int, std::vector<int> > connectivities;

        reader.ReadMesh(name, coordinates, connectivities);

        std::string corrected_name;
        bool found_mesh = false;
        if (mesh_info.Has(name))
        {
            corrected_name = name;
            found_mesh = true;
        }
        else
        {
            if (name.front() == '\"' && name.back() == '\"')
            {
                corrected_name = name.substr(1, name.size()-2);
                if (mesh_info.Has(corrected_name))
                    found_mesh = true;
            }
        }

        // get the entity name and type
        int entity_type, prop_id;
        std::string entity_name;
        bool do_guess_entity = true, do_guess_prop_id = true;

        if (found_mesh)
        {
            const Parameters& sub_params = mesh_info.GetValue(corrected_name);
            if (sub_params.Has("type"))
            {
                if (sub_params.GetValue("type").GetString() == "element")
                    entity_type = 1;
                else if (sub_params.GetValue("type").GetString() == "condition")
                    entity_type = 2;
                else
                    KRATOS_ERROR << "Invalid entity type " << sub_params.GetValue("type").GetString();
            }
            else
                entity_type = 0;

            if (sub_params.Has("name"))
            {
                entity_name = sub_params.GetValue("name").GetString();
                do_guess_entity = false;
            }

            if (sub_params.Has("prop_id"))
            {
                prop_id = sub_params.GetValue("prop_id").GetInt();
                do_guess_prop_id = false;
            }
        }

        if (do_guess_entity)
        {
            // try to guess the entity name from mesh name
            // TODO add more guess for element and condition

            if (name.find("Hexahedra3D27") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "DummyElement3D27N";  // to use this element, one needs to import the structural application. TODO move this element to kernel?
            }
            else if (name.find("Quadrilateral3D9") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "DummySurfaceCondition3D9N";
            }
            else
                KRATOS_ERROR << "Can't determine entity name for mesh " << name;
        }

        if (do_guess_prop_id)
        {
            // guess the prop_id from mesh name
            std::size_t pos = name.rfind('_');
            prop_id = std::stoi(name.substr(pos + 1));
        }

        if (echo_level > 0)
        {
            if (entity_type == 1)
                std::cout << "GiDPost2ModelPart: mesh " << name << " is assigned with element " << entity_name << " on prop_id " << prop_id << std::endl;
            else if (entity_type == 2)
                std::cout << "GiDPost2ModelPart: mesh " << name << " is assigned with condition " << entity_name << " on prop_id " << prop_id << std::endl;
        }

        // create nodes
        std::size_t last_node_id = r_model_part.GetLastNodeId();
        for (const auto& [id, coords] : coordinates)
        {
            r_model_part.CreateNewNode(id + last_node_id, coords[0], coords[1], coords[2]);
        }

        Properties::Pointer prop = r_model_part.pGetProperties(prop_id);

        if (entity_type == 1)
        {
            // create elements
            std::size_t last_elem_id = r_model_part.GetLastElementId();
            for (const auto& [id, conn] : connectivities)
            {
                std::vector<IndexType> new_conn;
                for (const int nid : conn)
                    new_conn.push_back(static_cast<IndexType>(nid + last_node_id));
                r_model_part.CreateNewElement(entity_name, id + last_elem_id, new_conn, prop);
            }
        }
        else if (entity_type == 2)
        {
            // create conditions
            std::size_t last_cond_id = r_model_part.GetLastConditionId();
            for (const auto& [id, conn] : connectivities)
            {
                std::vector<IndexType> new_conn;
                for (const int nid : conn)
                    new_conn.push_back(static_cast<IndexType>(nid + last_node_id));
                r_model_part.CreateNewCondition(entity_name, id + last_cond_id, new_conn, prop);
            }
        }
        else
            KRATOS_ERROR << "Entity type " << entity_type << " is invalid. Maybe you forgot to set the entity type?";

        // import nodal scalar results
        std::vector<std::string> nodal_scalar_result_names = reader.GetNodalScalarValuesName();

        const unsigned int buffer_size = r_model_part.GetBufferSize();
        if (echo_level > 0)
        {
            std::cout << "GiDPost2ModelPart: model_part buffer size = " << buffer_size << std::endl;
        }

        if (echo_level > 1)
        {
            if (nodal_scalar_result_names.size() > 0)
            {
                std::cout << "GiDPost2ModelPart: found nodal scalar values:";
                for (const auto& str : nodal_scalar_result_names)
                    std::cout << " " << str;
                std::cout << std::endl;
            }
            else
                std::cout << "GiDPost2ModelPart: found no nodal scalar values" << std::endl;
        }

        // import nodal vector results
        std::vector<std::string> nodal_vector_result_names = reader.GetNodalVectorValuesName();

        if (echo_level > 1)
        {
            if (nodal_vector_result_names.size() > 0)
            {
                std::cout << "GiDPost2ModelPart: found nodal vector values:";
                for (const auto& str : nodal_vector_result_names)
                    std::cout << " " << str;
                std::cout << std::endl;
            }
            else
                std::cout << "GiDPost2ModelPart: found no nodal vector values" << std::endl;
        }

        if (nodal_vector_result_names.size() > 0)
        {
            const int vector_size = 4; // for now we support reading vector values of size 3 only TODO generalize the vector dimension.

            for (const auto& result_name : nodal_vector_result_names)
            {
                // read the nodal values
                std::vector<double> step_list;
                std::map<std::size_t, std::vector<std::vector<double> > > step_values;

                reader.ReadNodalVectorValues(result_name, step_list, step_values, vector_size);

                // get the variable from kernel
                std::string var_name = GiDPostReader::StripQuote(result_name);
                if (!KratosComponents<Variable<array_1d<double, 3> > >::Has(var_name))
                    KRATOS_ERROR << "Variable " << var_name << " is not registerred to the kernel";
                const Variable<array_1d<double, 3> >& rVariable = KratosComponents<Variable<array_1d<double, 3> > >::Get(var_name);

                for (const auto& [id, values] : step_values)
                {
                    auto& rNode = r_model_part.Nodes()[id + last_node_id];

                    unsigned int cnt = 0;
                    for (auto it = values.rbegin(); it != values.rend(); ++it)
                    {
                        if (cnt >= buffer_size)
                            break; // only read until buffer size
                        array_1d<double, 3> v;
                        v[0] = (*it)[0];
                        v[1] = (*it)[1];
                        v[2] = (*it)[2];

                        noalias(rNode.GetSolutionStepValue(rVariable, cnt)) = v;

                        if (echo_level > 2)
                        {
                            std::cout << "GiDPost2ModelPart: node " << rNode.Id() << " is assigned " << rVariable.Name()
                                      << " with value " << v << " at position " << cnt
                                      << std::endl;
                        }

                        ++cnt;
                    }
                }

                if (echo_level > 0)
                {
                    std::cout << "GiDPost2ModelPart: read nodal vector values " << rVariable.Name() << " completed" << std::endl;
                }
            }
        }
    }
}

}// namespace Kratos.
