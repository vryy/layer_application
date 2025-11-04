//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 22 Jul 2019 $
//   Revision:            $Revision: 1.0 $
//
//
// Project includes
#include "includes/legacy_structural_app_vars.h"
#include "layer_application_variables.h"
#if defined(PROFILING_LEVEL)
#include "utilities/openmp_utils.h"
#endif
#include "utilities/math_utils.h"
#include "utilities/string_utils.h"
#include "utilities/triangulation_utils.h"
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

void ModelPartUtilities::ExportEdgeInformation(std::ostream& rOStream, const ModelPart::ElementsContainerType& rpElements,
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

void ModelPartUtilities::ExportEdgeInformationToGiD(std::ostream& rOStream, const ModelPart::ElementsContainerType& rpElements)
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

int ModelPartUtilities::GetIntegrationOrder(const GeometryType& rGeometry, const unsigned int npoints)
{
    unsigned int number_of_integration_methods = static_cast<unsigned int>(GeometryData::IntegrationMethod::NumberOfIntegrationMethods);
    for (unsigned int i = 0; i < number_of_integration_methods; ++i)
    {
        const GeometryData::IntegrationMethod ThisIntegrationMethod =
            static_cast<GeometryData::IntegrationMethod>(i+1);

        const GeometryType::IntegrationPointsArrayType& integration_points =
            rGeometry.IntegrationPoints( ThisIntegrationMethod );

        if (integration_points.size() == npoints)
            return i + 1;
    }

    KRATOS_ERROR << "Can't find suitable integration order for " << npoints
                 << " integration points of geometry " << rGeometry.GetGeometryType();
}

// TODO to test
void ModelPartUtilities::ClearModelPart(ModelPart& r_model_part)
{
    // remove all elements
    std::vector<IndexType> all_elements;
    for (auto it = r_model_part.ElementsBegin(); it != r_model_part.ElementsEnd(); ++it)
        all_elements.push_back(it->Id());
    for (const auto id : all_elements)
        r_model_part.RemoveElement(id);

    // remove all conditions
    std::vector<IndexType> all_conditions;
    for (auto it = r_model_part.ConditionsBegin(); it != r_model_part.ConditionsEnd(); ++it)
        all_conditions.push_back(it->Id());
    for (const auto id : all_conditions)
        r_model_part.RemoveCondition(id);

    // remove all nodes
    std::vector<IndexType> all_nodes;
    for (auto it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it)
        all_nodes.push_back(it->Id());
    for (const auto id : all_nodes)
        r_model_part.RemoveNode(id);

    // remove all properties
    std::vector<IndexType> all_properties;
    for (auto it = r_model_part.PropertiesBegin(); it != r_model_part.PropertiesEnd(); ++it)
        all_properties.push_back(it->Id());
    for (const auto id : all_properties)
        r_model_part.RemoveProperties(id);

    // remove all constraints
    std::vector<IndexType> all_constraints;
    for (auto it = r_model_part.MasterSlaveConstraintsBegin(); it != r_model_part.MasterSlaveConstraintsEnd(); ++it)
        all_constraints.push_back(it->Id());
    for (const auto id : all_constraints)
        r_model_part.RemoveMasterSlaveConstraint(id);
}

void ModelPartUtilities::CreateSubModelPartFromElements(ModelPart& r_model_part, const std::string& Name, const std::vector<std::size_t>& element_ids)
{
    ModelPart& r_sub_model_part = dynamic_cast<ModelPart&>(r_model_part.CreateSubModelPart(Name));

    std::set<ModelPart::IndexType> node_set;

    auto& rElements = r_model_part.Elements();
    for (auto elem_id : element_ids)
    {
        auto elem = ModelPart::ElementType::Pointer(rElements(elem_id));
        r_sub_model_part.AddElement(elem);

        const auto& rGeometry = elem->GetGeometry();
        for (std::size_t i = 0; i < rGeometry.size(); ++i)
            node_set.insert(rGeometry[i].Id());
    }

    auto& rNodes = r_model_part.Nodes();
    for (auto node_id : node_set)
    {
        auto node = ModelPart::NodeType::Pointer(rNodes(node_id));
        r_sub_model_part.AddNode(node);
    }
}

void ModelPartUtilities::GiDPostBin2ModelPart(const std::string& fileName, ModelPart& r_model_part,
        const Parameters& mesh_info, VariablesList<>* pElementalVariablesList)
{
    GiDPostBinaryReader reader(fileName);
    GiDPost2ModelPart(reader, r_model_part, mesh_info, pElementalVariablesList);
    r_model_part.Name() = reader.GetPrefix();               // set the name
    r_model_part.GetProcessInfo()[TIME] = reader.GetTime(); // set the time
}

void ModelPartUtilities::GiDPost2ModelPart(GiDPostReader& reader, ModelPart& r_model_part,
        const Parameters& mesh_info, VariablesList<>* pElementalVariablesList)
{
    typedef ModelPart::IndexType IndexType;

    std::vector<std::string> mesh_names = reader.GetMeshesName();

    int echo_level = mesh_info.Has("echo_level") ? mesh_info.GetValue("echo_level").GetInt() : 0;

    // we iterate through all mesh to collect all nodal variables.
    // This ensures all new nodes are created consistently with all
    // nodal variables properly set
    {
        // here we have to get all nodal result to add the solution step variables to node
        std::vector<std::string> nodal_scalar_result_names = reader.GetNodalScalarValuesName();
        std::vector<std::string> nodal_vector_result_names = reader.GetNodalVectorValuesName();

        for (const auto& result_name : nodal_scalar_result_names)
        {
            // get the variable from kernel
            std::string var_name = StringUtils::StripQuote(result_name, '\"');
            if (!KratosComponents<Variable<double> >::Has(var_name))
                KRATOS_ERROR << "Variable " << var_name << " is not registerred to the kernel";
            const Variable<double>& rVariable = KratosComponents<Variable<double> >::Get(var_name);

            if (!r_model_part.HasNodalSolutionStepVariable(rVariable))
                r_model_part.AddNodalSolutionStepVariable(rVariable);

            if (echo_level > 1)
                std::cout << "GiDPost2ModelPart: found nodal scalar value " << result_name << std::endl;
        }

        for (const auto& result_name : nodal_vector_result_names)
        {
            // get the variable from kernel
            std::string var_name = StringUtils::StripQuote(result_name, '\"');
            if (!KratosComponents<Variable<array_1d<double, 3> > >::Has(var_name))
                KRATOS_ERROR << "Variable " << var_name << " is not registerred to the kernel";
            const Variable<array_1d<double, 3> >& rVariable = KratosComponents<Variable<array_1d<double, 3> > >::Get(var_name);

            if (!r_model_part.HasNodalSolutionStepVariable(rVariable))
                r_model_part.AddNodalSolutionStepVariable(rVariable);

            if (echo_level > 1)
                std::cout << "GiDPost2ModelPart: found nodal vector value " << result_name << std::endl;
        }
    }

    // we append the node, element and condition from the corresponding last index of the input model_part
    IndexType last_node_id = r_model_part.GetLastNodeId();
    IndexType last_elem_id = r_model_part.GetLastElementId();
    IndexType last_cond_id = r_model_part.GetLastConditionId();

    /* create nodes */

    for (const auto& name : mesh_names)
    {
        std::map<int, std::vector<double> > coordinates;

#if PROFILING_LEVEL > 0
        double time1 = OpenMPUtils::GetCurrentTime();
#endif

        reader.ReadMesh(name, coordinates);

#if PROFILING_LEVEL > 0
        double time2 = OpenMPUtils::GetCurrentTime();
#endif

        SizeType nnodes = 0;
        for (const auto& [id, coords] : coordinates)
        {
            // because the node can be duplicated. Note that each mesh
            // can have more than nodes that it required, just the element
            // and condition are unique.
            if (!r_model_part.HasNode(id + last_node_id))
            {
                r_model_part.CreateNewNode(id + last_node_id, coords[0], coords[1], coords[2]);
                ++nnodes;

                if (echo_level > 2)
                {
                    std::cout << "GiDPost2ModelPart: node " << id + last_node_id << " is created"
                              << std::endl;
                }
            }
        }

#if PROFILING_LEVEL > 0
        double time3 = OpenMPUtils::GetCurrentTime();
#endif

        if (echo_level > 1)
            std::cout << "GiDPost2ModelPart: " << nnodes << " nodes are created on mesh " << name
#if PROFILING_LEVEL > 0
                      << ", read coordinates = " << (time2 - time1) << "s"
                      << ", create nodes = " << (time3 - time2) << "s"
#endif
                      << std::endl;
    } // end creating nodes

    /* import nodal scalar results */

    const unsigned int buffer_size = r_model_part.GetBufferSize();
    if (echo_level > 0)
    {
        std::cout << "GiDPost2ModelPart: model_part buffer size = " << buffer_size << std::endl;
    }

    std::vector<std::string> nodal_scalar_result_names = reader.GetNodalScalarValuesName();

    for (const auto& result_name : nodal_scalar_result_names)
    {
        // read the nodal values
        std::vector<double> step_list;
        std::map<std::size_t, std::vector<double> > step_values;

#if PROFILING_LEVEL > 0
        double time1 = OpenMPUtils::GetCurrentTime();
#endif

        reader.ReadNodalScalarValues(result_name, step_list, step_values);

#if PROFILING_LEVEL > 0
        double time2 = OpenMPUtils::GetCurrentTime();
#endif

        // get the variable from kernel
        std::string var_name = StringUtils::StripQuote(result_name, '\"');
        const Variable<double>& rVariable = KratosComponents<Variable<double> >::Get(var_name);

        for (const auto& [id, values] : step_values)
        {
            auto it = r_model_part.Nodes().find(id + last_node_id);
            if (it == r_model_part.Nodes().end())
                continue;

            auto& rNode = *it;

            unsigned int cnt = 0;
            for (auto it = values.rbegin(); it != values.rend(); ++it)
            {
                if (cnt >= buffer_size)
                    break; // only read until buffer size

                rNode.GetSolutionStepValue(rVariable, cnt) = *it;

                if (echo_level > 2)
                {
                    std::cout << "GiDPost2ModelPart: node " << rNode.Id() << " is assigned " << rVariable.Name()
                              << " with value " << *it << " at position " << cnt
                              << std::endl;
                }

                ++cnt;
            }
        }

#if PROFILING_LEVEL > 0
        double time3 = OpenMPUtils::GetCurrentTime();
#endif

        if (echo_level > 1)
        {
            std::cout << "GiDPost2ModelPart: read nodal scalar values " << rVariable.Name() << " completed"
#if PROFILING_LEVEL > 0
                      << ", read results: " << (time2 - time1) << "s"
                      << ", transfer results: " << (time3 - time2) << "s"
#endif
                      << std::endl;
        }
    } // end importing nodal scalar results

    /* import nodal vector results */

    std::vector<std::string> nodal_vector_result_names = reader.GetNodalVectorValuesName();

    for (const auto& result_name : nodal_vector_result_names)
    {
        const int vector_size = 4; // for now we support reading vector values of size 3 only
        // TODO generalize the vector dimension., maybe through json parameters?

        // read the nodal values
        std::vector<double> step_list;
        std::map<std::size_t, std::vector<std::vector<double> > > step_values;

#if PROFILING_LEVEL > 0
        double time1 = OpenMPUtils::GetCurrentTime();
#endif

        reader.ReadNodalVectorValues(result_name, step_list, step_values, vector_size);

#if PROFILING_LEVEL > 0
        double time2 = OpenMPUtils::GetCurrentTime();
#endif

        // get the variable from kernel
        std::string var_name = StringUtils::StripQuote(result_name, '\"');
        const Variable<array_1d<double, 3> >& rVariable = KratosComponents<Variable<array_1d<double, 3> > >::Get(var_name);

        for (const auto& [id, values] : step_values)
        {
            auto it = r_model_part.Nodes().find(id + last_node_id);
            if (it == r_model_part.Nodes().end())
                continue;

            auto& rNode = *it;

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

#if PROFILING_LEVEL > 0
        double time3 = OpenMPUtils::GetCurrentTime();
#endif

        if (echo_level > 1)
        {
            std::cout << "GiDPost2ModelPart: read nodal vector values " << rVariable.Name() << " completed"
#if PROFILING_LEVEL > 0
                      << ", read results: " << (time2 - time1) << "s"
                      << ", transfer results: " << (time3 - time2) << "s"
#endif
                      << std::endl;
        }
    } // end importing nodal vector results

    /* create elements and conditions */

    for (const auto& name : mesh_names)
    {
        std::map<int, std::vector<int> > connectivities;

        reader.ReadMesh(name, connectivities);

        int dim;
        std::string gid_elem_type;
        reader.GetMeshInfo(name, dim, gid_elem_type);

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
        std::string layer_name = "", entity_name;
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

            if (sub_params.Has("layer_name"))
            {
                layer_name = sub_params.GetValue("layer_name").GetString();
            }
        }

        if (do_guess_entity)
        {
            // try to guess the entity name from mesh name
            // TODO add more guess for element and condition

            if (name.find("Hexahedra3D27") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement3D27N";
            }
            else if (name.find("Hexahedra3D20") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement3D20N";
            }
            else if (name.find("Hexahedra3D8") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement3D8N";
            }
            else if (name.find("Tetrahedra3D10") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement3D10N";
            }
            else if (name.find("Tetrahedra3D4") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement3D4N";
            }
            else if (name.find("Quadrilateral2D9") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement2D9N";
            }
            else if (name.find("Quadrilateral2D8") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement2D8N";
            }
            else if (name.find("Quadrilateral2D4") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement2D4N";
            }
            else if (name.find("Triangle2D6") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement2D6N";
            }
            else if (name.find("Triangle2D3") != std::string::npos)
            {
                entity_type = 1;
                entity_name = "PostElement2D3N";
            }
            else if (name.find("Quadrilateral3D9") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostSurfaceCondition3D9N";
            }
            else if (name.find("Quadrilateral3D8") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostSurfaceCondition3D8N";
            }
            else if (name.find("Quadrilateral3D4") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostSurfaceCondition3D4N";
            }
            else if (name.find("Triangle3D6") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostSurfaceCondition3D6N";
            }
            else if (name.find("Triangle3D3") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostSurfaceCondition3D3N";
            }
            else if (name.find("Line2D2") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostLineCondition2D2N";
            }
            else if (name.find("Line2D3") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostLineCondition2D3N";
            }
            else if (name.find("Line3D2") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostLineCondition3D2N";
            }
            else if (name.find("Line3D3") != std::string::npos)
            {
                entity_type = 2;
                entity_name = "PostLineCondition3D3N";
            }
            else
                KRATOS_ERROR << "Can't determine entity name for mesh " << name;
        }

        if (do_guess_prop_id)
        {
            // guess the prop_id from mesh name
            ExtractPropertiesId(name, prop_id, layer_name);
            if (prop_id < 0)
                // KRATOS_ERROR << "Failed to extract Properties Id from mesh " << name;
                prop_id = 1;

            if (echo_level > 1)
            {
                std::cout << "GiDPost2ModelPart: Properties Id " << prop_id << " is predicted for mesh " << name
                          << std::endl;
            }
        }

        // do some checks
        if (dim == 2)
        {
            if ((gid_elem_type == "Quadrilateral" || gid_elem_type == "Triangle")
                && entity_type != 1
            )
                KRATOS_ERROR << "Mesh is " << gid_elem_type << " but the entity is not element. Is there something wrong?";

            if ((gid_elem_type == "Line" || gid_elem_type == "Linear")
                && entity_type != 2
            )
                KRATOS_ERROR << "Mesh is " << gid_elem_type << " but the entity is not condition. Is there something wrong?";
        }
        if (dim == 3)
        {
            if ((gid_elem_type == "Hexahedra" || gid_elem_type == "Tetrathedra")
                && entity_type != 1
            )
                KRATOS_ERROR << "Mesh is " << gid_elem_type << " but the entity is not element. Is there something wrong?";

            if ((gid_elem_type == "Quadrilateral" || gid_elem_type == "Triangle")
                && entity_type != 2
            )
                KRATOS_ERROR << "Mesh is " << gid_elem_type << " but the entity is not condition. Is there something wrong?";
        }

        if (echo_level > 0)
        {
            if (entity_type == 1)
                std::cout << "GiDPost2ModelPart: mesh " << name << " is assigned with element " << entity_name << " on prop_id " << prop_id << std::endl;
            else if (entity_type == 2)
                std::cout << "GiDPost2ModelPart: mesh " << name << " is assigned with condition " << entity_name << " on prop_id " << prop_id << std::endl;
        }

        Properties::Pointer prop = r_model_part.pGetProperties(prop_id);
        prop->SetValue(LAYER_NAME, layer_name);

        if (entity_type == 1)
        {
            // create elements
            for (const auto& [id, conn] : connectivities)
            {
                std::vector<IndexType> new_conn;
                for (const int nid : conn)
                    new_conn.push_back(static_cast<IndexType>(nid + last_node_id));
                r_model_part.CreateNewElement(entity_name, id + last_elem_id, new_conn, prop)
                ->Initialize(r_model_part.GetProcessInfo());
            }
            if (echo_level > 1)
                std::cout << "GiDPost2ModelPart: " << connectivities.size() << " elements are created on mesh " << name << std::endl;
        }
        else if (entity_type == 2)
        {
            // create conditions
            for (const auto& [id, conn] : connectivities)
            {
                std::vector<IndexType> new_conn;
                for (const int nid : conn)
                    new_conn.push_back(static_cast<IndexType>(nid + last_node_id));
                r_model_part.CreateNewCondition(entity_name, id + last_cond_id, new_conn, prop)
                ->Initialize(r_model_part.GetProcessInfo());
            }
            if (echo_level > 1)
                std::cout << "GiDPost2ModelPart: " << connectivities.size() << " conditions are created on mesh " << name << std::endl;
        }
        else
            KRATOS_ERROR << "Entity type " << entity_type << " is invalid. Maybe you forgot to set the entity type?";

        /* import Gauss point scalar results */

        std::vector<std::pair<std::string, std::string> > gp_scalar_result_names = reader.GetGaussPointScalarValuesName();

        for (const auto& [result_name, gp_name] : gp_scalar_result_names)
        {
            // get the variable from kernel
            std::string var_name = StringUtils::StripQuote(result_name, '\"');
            if (!KratosComponents<Variable<double> >::Has(var_name))
            {
                // KratosComponents<Variable<double> >::Add(var_name) // TODO to create new variable and add to the kernel
                std::cout << "WARNING!!!Variable " << var_name << " is not registerred to the kernel. The result associated with it is skipped.";
                continue;
            }
            const Variable<double>& rVariable = KratosComponents<Variable<double> >::Get(var_name);

            // get the Gauss point record info
            int npoints;
            std::string gp_elem_type;
            std::string gp_coordinates_type;
            reader.GetGaussPointRecordInfo(gp_name, npoints, gp_elem_type, gp_coordinates_type);
            reader.ReadGaussPointRecord(gp_name);

            std::vector<array_1d<double, 3> > gp_coordinates;
            reader.GetGaussPointRecordCoordinates(gp_name, gp_coordinates);

            // results (with name) can be on different mesh. It is essential to check if the element type
            // of the Gauss point results match with the mesh element type
            if (gid_elem_type.compare(gp_elem_type) != 0)
                continue;

            if (echo_level > 1)
                std::cout << "GiDPost2ModelPart: found Gauss point scalar value " << result_name
                          << ", name " << gp_name
                          << std::endl;

            if (pElementalVariablesList != nullptr)
            {
                pElementalVariablesList->Add(rVariable);
            }

            std::vector<double> step_list;
            std::map<std::size_t, std::vector<std::vector<double> > > step_values;
            reader.ReadGaussPointScalarValues(result_name, gp_name, step_list, step_values);

            for (const auto& [id, values] : step_values)
            {
                if (entity_type == 1)
                {
                    auto it = r_model_part.Elements().find(id + last_elem_id);
                    if (it == r_model_part.ElementsEnd())
                        continue;

                    auto& rElement = *it;

                    for (auto it2 = values.rbegin(); it2 != values.rend(); ++it2)
                    {
                        // do a size check
                        if (it2->size() != npoints)
                            KRATOS_ERROR << "The number of values does not match number of integration points in the Gauss point record. Is there something wrong?";

                        if (gp_coordinates_type == "Given")
                        {
                            // set the integration point coordinates. The PostElement will need it to calculate the global coordinates of the integration point
                            rElement.SetValuesOnIntegrationPoints(INTEGRATION_POINT_LOCAL, gp_coordinates, r_model_part.GetProcessInfo());
                        }
                        else if (gp_coordinates_type == "Internal")
                        {
                            // guess the integration order based on number of integration points
                            int integration_order = GetIntegrationOrder(rElement.GetGeometry(), npoints);
                            if (rElement.Has(INTEGRATION_ORDER))
                            {
                                if (rElement.GetValue(INTEGRATION_ORDER) != integration_order)
                                    KRATOS_ERROR << "The existing INTEGRATION_ORDER of element " << rElement.Id()
                                                 << " conflicts with the number of integration points on Gp record."
                                                 << " Something is inconsistent.";
                            }
                            else
                                rElement.SetValue(INTEGRATION_ORDER, integration_order);
                        }
                        else
                            KRATOS_ERROR << "Unknown coordinates type " << gp_coordinates_type;

                        rElement.SetValuesOnIntegrationPoints(rVariable, *it2, r_model_part.GetProcessInfo());

                        if (echo_level > 2)
                        {
                            std::cout << "GiDPost2ModelPart: element " << rElement.Id() << " is assigned " << rVariable.Name()
                                      << " with values";
                            for (std::size_t i = 0; i < it2->size(); ++i)
                                std::cout << " " << (*it2)[i];
                            std::cout << std::endl;
                        }
                    }
                }
                else if (entity_type == 2)
                {
                    auto it = r_model_part.Conditions().find(id + last_cond_id);
                    if (it == r_model_part.ConditionsEnd())
                        continue;

                    auto& rCondition = *it;

                    for (auto it2 = values.rbegin(); it2 != values.rend(); ++it2)
                    {
                        // do a size check
                        if (it2->size() != npoints)
                            KRATOS_ERROR << "The number of values does not match number of integration points in the Gauss point record. Is there something wrong?";

                        if (gp_coordinates_type == "Given")
                        {
                            // set the integration point coordinates. The PostCondition will need it to calculate the global coordinates of the integration point
                            rCondition.SetValuesOnIntegrationPoints(INTEGRATION_POINT_LOCAL, gp_coordinates, r_model_part.GetProcessInfo());
                        }
                        else if (gp_coordinates_type == "Internal")
                        {
                            // guess the integration order based on number of integration points
                            int integration_order = GetIntegrationOrder(rCondition.GetGeometry(), npoints);
                            if (rCondition.Has(INTEGRATION_ORDER))
                            {
                                if (rCondition.GetValue(INTEGRATION_ORDER) != integration_order)
                                    KRATOS_ERROR << "The existing INTEGRATION_ORDER of element " << rCondition.Id()
                                                 << " conflicts with the number of integration points on Gp record."
                                                 << " Something is inconsistent.";
                            }
                            else
                                rCondition.SetValue(INTEGRATION_ORDER, integration_order);
                        }
                        else
                            KRATOS_ERROR << "Unknown coordinates type " << gp_coordinates_type;

                        rCondition.SetValuesOnIntegrationPoints(rVariable, *it2, r_model_part.GetProcessInfo());

                        if (echo_level > 2)
                        {
                            std::cout << "GiDPost2ModelPart: condition " << rCondition.Id() << " is assigned " << rVariable.Name()
                                      << " with values";
                            for (std::size_t i = 0; i < it2->size(); ++i)
                                std::cout << " " << (*it2)[i];
                            std::cout << std::endl;
                        }
                    }
                }
            }

            if (echo_level > 0)
            {
                std::cout << "GiDPost2ModelPart: read Gauss point scalar values " << rVariable.Name() << " completed" << std::endl;
            }
        }

        /* import Gauss point vector results */

        std::vector<std::pair<std::string, std::string> > gp_vector_result_names = reader.GetGaussPointVectorValuesName();

        for (const auto& [result_name, gp_name] : gp_vector_result_names)
        {
            // get the variable from kernel
            std::string var_name = StringUtils::StripQuote(result_name, '\"');
            if (!KratosComponents<Variable<Vector> >::Has(var_name))
            {
                // KratosComponents<Variable<Vector> >::Add(var_name) // TODO to create new variable and add to the kernel
                std::cout << "WARNING!!!Variable " << var_name << " is not registerred to the kernel. The result associated with it is skipped.";
                continue;
            }
            const Variable<Vector>& rVariable = KratosComponents<Variable<Vector> >::Get(var_name);

            // get the Gauss point record info
            int npoints;
            std::string gp_elem_type;
            std::string gp_coordinates_type;
            reader.GetGaussPointRecordInfo(gp_name, npoints, gp_elem_type, gp_coordinates_type);
            reader.ReadGaussPointRecord(gp_name);

            std::vector<array_1d<double, 3> > gp_coordinates;
            reader.GetGaussPointRecordCoordinates(gp_name, gp_coordinates);

            // results (with name) can be on different mesh. It is essential to check if the element type
            // of the Gauss point results match with the mesh element type
            if (gid_elem_type.compare(gp_elem_type) != 0)
                continue;

            if (echo_level > 1)
                std::cout << "GiDPost2ModelPart: found Gauss point vector value " << result_name
                          << ", name " << gp_name
                          << std::endl;

            if (pElementalVariablesList != nullptr)
            {
                pElementalVariablesList->Add(rVariable);
            }

            // TODO to read the value and assign to element integration points
        }

        /* import Gauss point matrix results */

        std::vector<std::pair<std::string, std::string> > gp_matrix_result_names = reader.GetGaussPointMatrixValuesName();

        for (const auto& [result_name, gp_name] : gp_matrix_result_names)
        {
            // get the variable from kernel
            std::string var_name = StringUtils::StripQuote(result_name, '\"');

            if (KratosComponents<Variable<Vector> >::Has(var_name))
            // some matrix results are written as vector. Here we try to read the matrix results as vector first
            {
                const Variable<Vector>& rVariable = KratosComponents<Variable<Vector> >::Get(var_name);

                std::cout << "Matrix result " << var_name << " is detected to save to vector varable" << std::endl;

                // get the Gauss point record info
                int npoints;
                std::string gp_elem_type;
                std::string gp_coordinates_type;
                reader.GetGaussPointRecordInfo(gp_name, npoints, gp_elem_type, gp_coordinates_type);
                reader.ReadGaussPointRecord(gp_name);

                std::vector<array_1d<double, 3> > gp_coordinates;
                reader.GetGaussPointRecordCoordinates(gp_name, gp_coordinates);

                // results (with name) can be on different mesh. It is essential to check if the element type
                // of the Gauss point results match with the mesh element type
                if (gid_elem_type.compare(gp_elem_type) != 0)
                    continue;

                if (echo_level > 1)
                    std::cout << "GiDPost2ModelPart: found Gauss point matrix value " << result_name
                              << ", name " << gp_name
                              << std::endl;

                if (pElementalVariablesList != nullptr)
                {
                    pElementalVariablesList->Add(rVariable);
                }

                std::vector<double> step_list;
                std::map<std::size_t, std::vector<std::vector<std::vector<double> > > > step_values;
                reader.ReadGaussPointMatrixValues(result_name, gp_name, step_list, step_values);

                for (const auto& [id, values] : step_values)
                {
                    if (entity_type == 1)
                    {
                        auto it = r_model_part.Elements().find(id + last_elem_id);
                        if (it == r_model_part.ElementsEnd())
                            continue;

                        auto& rElement = *it;

                        for (auto it2 = values.rbegin(); it2 != values.rend(); ++it2)
                        {
                            // do a size check
                            if (it2->size() != npoints)
                                KRATOS_ERROR << "The number of values does not match number of integration points in the Gauss point record. Is there something wrong?";

                            if (gp_coordinates_type == "Given")
                            {
                                // set the integration point coordinates. The PostElement will need it to calculate the global coordinates of the integration point
                                rElement.SetValuesOnIntegrationPoints(INTEGRATION_POINT_LOCAL, gp_coordinates, r_model_part.GetProcessInfo());
                            }
                            else if (gp_coordinates_type == "Internal")
                            {
                                // guess the integration order based on number of integration points
                                int integration_order = GetIntegrationOrder(rElement.GetGeometry(), npoints);
                                if (rElement.Has(INTEGRATION_ORDER))
                                {
                                    if (rElement.GetValue(INTEGRATION_ORDER) != integration_order)
                                        KRATOS_ERROR << "The existing INTEGRATION_ORDER of element " << rElement.Id()
                                                     << " conflicts with the number of integration points on Gp record."
                                                     << " Something is inconsistent.";
                                }
                                else
                                    rElement.SetValue(INTEGRATION_ORDER, integration_order);
                            }
                            else
                                KRATOS_ERROR << "Unknown coordinates type " << gp_coordinates_type;

                            std::vector<Vector> temp(it2->size());
                            auto it3 = it2->begin();
                            for (std::size_t i = 0; i < it2->size(); ++i, ++it3)
                            {
                                temp[i].resize(it3->size());
                                std::copy(it3->begin(), it3->end(), temp[i].begin());
                            }
                            rElement.SetValuesOnIntegrationPoints(rVariable, temp, r_model_part.GetProcessInfo());

                            if (echo_level > 2)
                            {
                                std::cout << "GiDPost2ModelPart: element " << rElement.Id() << " is assigned " << rVariable.Name()
                                          << " with values";
                                for (std::size_t i = 0; i < it2->size(); ++i)
                                {
                                    for (std::size_t j = 0; j < (*it2)[i].size(); ++j)
                                        std::cout << " " << (*it2)[i][j];
                                    std::cout << ";";
                                }
                                std::cout << std::endl;
                            }
                        }
                    }
                    else if (entity_type == 2)
                    {
                        auto it = r_model_part.Conditions().find(id + last_cond_id);
                        if (it == r_model_part.ConditionsEnd())
                            continue;

                        auto& rCondition = *it;

                        for (auto it2 = values.rbegin(); it2 != values.rend(); ++it2)
                        {
                            // do a size check
                            if (it2->size() != npoints)
                                KRATOS_ERROR << "The number of values does not match number of integration points in the Gauss point record. Is there something wrong?";

                            if (gp_coordinates_type == "Given")
                            {
                                // set the integration point coordinates. The PostCondition will need it to calculate the global coordinates of the integration point
                                rCondition.SetValuesOnIntegrationPoints(INTEGRATION_POINT_LOCAL, gp_coordinates, r_model_part.GetProcessInfo());
                            }
                            else if (gp_coordinates_type == "Internal")
                            {
                                // guess the integration order based on number of integration points
                                int integration_order = GetIntegrationOrder(rCondition.GetGeometry(), npoints);
                                if (rCondition.Has(INTEGRATION_ORDER))
                                {
                                    if (rCondition.GetValue(INTEGRATION_ORDER) != integration_order)
                                        KRATOS_ERROR << "The existing INTEGRATION_ORDER of element " << rCondition.Id()
                                                     << " conflicts with the number of integration points on Gp record."
                                                     << " Something is inconsistent.";
                                }
                                else
                                    rCondition.SetValue(INTEGRATION_ORDER, integration_order);
                            }
                            else
                                KRATOS_ERROR << "Unknown coordinates type " << gp_coordinates_type;

                            std::vector<Vector> temp(it2->size());
                            auto it3 = it2->begin();
                            for (std::size_t i = 0; i < it2->size(); ++i, ++it3)
                            {
                                temp[i].resize(it3->size());
                                std::copy(it3->begin(), it3->end(), temp[i].begin());
                            }
                            rCondition.SetValuesOnIntegrationPoints(rVariable, temp, r_model_part.GetProcessInfo());

                            if (echo_level > 2)
                            {
                                std::cout << "GiDPost2ModelPart: condition " << rCondition.Id() << " is assigned " << rVariable.Name()
                                          << " with values";
                                for (std::size_t i = 0; i < it2->size(); ++i)
                                {
                                    for (std::size_t j = 0; j < (*it2)[i].size(); ++j)
                                        std::cout << " " << (*it2)[i][j];
                                    std::cout << ";";
                                }
                                std::cout << std::endl;
                            }
                        }
                    }
                }

                if (echo_level > 0)
                {
                    std::cout << "GiDPost2ModelPart: read Gauss point matrix values as vector variable " << rVariable.Name() << " completed" << std::endl;
                }
            }
            // else if (KratosComponents<Variable<Matrix> >::Has(var_name))
            // {
            //     // TODO
            // }
            else
            {
                std::cout << "WARNING!!!Variable " << var_name << " is not registerred to the kernel. The result associated with it is skipped.";
                continue;
            }
        }
    } // end creating elements and conditions
}

void ModelPartUtilities::ExtractPropertiesId(const std::string& name, int& prop_id, std::string& layer_name)
{
    // strip the quote if needed
    std::string name_correct = StringUtils::StripQuote(name, '\"');

    // split to words
    std::vector<std::string> words = StringUtils::Split(name_correct, '_');

    // // look back from last and check which one is an integer
    // for (auto it = words.rbegin(); it != words.rend(); ++it)
    // {
    //     if (StringUtils::IsValidInteger(*it))
    //         return std::atoi(it->c_str());
    // }

    // extract all the integers
    std::vector<int> numbers;
    std::map<int, std::size_t> map_number_to_index;
    std::size_t cnt = 0;
    for (auto it = words.begin(); it != words.end(); ++it)
    {
        if (StringUtils::IsValidInteger(*it))
        {
            int i = std::atoi(it->c_str());
            map_number_to_index[i] = cnt;
            numbers.push_back(i);
        }
        ++cnt;
    }

    if (numbers.size() > 0)
    {
        prop_id = numbers[0]; // the first number is chosen as Properties Id
        std::size_t i = map_number_to_index[prop_id];
        layer_name = "";
        for (std::size_t j = i+1; j < words.size(); ++j)
        {
            if (j == words.size() - 1)
                layer_name += words[j];
            else
                layer_name += words[j] + "_";
        }
        return;
    }

    prop_id = -1;
}

void ModelPartUtilities::XY2ModelPart(const std::vector<array_1d<CoordinateType, 3> >& rPoints, ModelPart& r_model_part,
            const IndexType prop_id, const array_1d<double, 3>& rDirector, const std::string& condition_name)
{
    const IndexType last_node_id = r_model_part.GetLastNodeId();
    const IndexType last_cond_id = r_model_part.GetLastConditionId();

    array_1d<double, 3> xDir, yDir;
    MathUtils<double>::ComputeLocalCoordinateSystem(rDirector, xDir, yDir);

    // add node to the model_part
    // and also create list of point projection to compute connectivity
    std::vector<double> xylist;
    IndexType new_node_id = last_node_id;
    for (const auto& p : rPoints)
    {
        r_model_part.CreateNewNode(++new_node_id, p[0], p[1], p[2]);

        xylist.push_back(inner_prod(p, xDir));
        xylist.push_back(inner_prod(p, yDir));
    }

    // triangulation
    std::vector<std::vector<IndexType> > connectivities;
    TriangulationUtils tri_utils;
    tri_utils.ComputeDelaunayTriangulation(xylist, connectivities);

    IndexType new_cond_id = last_cond_id;
    auto prop = r_model_part.pGetProperties(prop_id);
    for (const auto& c : connectivities)
    {
        std::vector<IndexType> conn = c;
        for (auto& i : conn)
            i += last_node_id + 1;

        r_model_part.CreateNewCondition(condition_name, ++new_cond_id, conn, prop);
    }
}

} // namespace Kratos.
