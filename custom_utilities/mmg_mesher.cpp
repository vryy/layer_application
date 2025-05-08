#include "mmg_mesher.h"
#include "layer_application/layer_application_variables.h"

namespace Kratos
{

template<int TDim>
MMGMesher<TDim>::MMGMesher()
{
    mmmgMesh = NULL;
    mmmgMet = NULL;
    mmmgLs = NULL;

    muse_level_set = false;
    muse_metric = true;

    MMG<TDim>::Init_mesh_met( mmmgMesh, mmmgMet );
}

template<int TDim>
MMGMesher<TDim>::MMGMesher(const bool& use_level_set, const bool& use_metric)
{
    mmmgMesh = NULL;
    mmmgMet = NULL;
    mmmgLs = NULL;

    muse_level_set = use_level_set;
    muse_metric = use_metric;

    if (!muse_level_set)
    {
        muse_metric = true; // if level set mode is not active, the metric is required
        MMG<TDim>::Init_mesh_met( mmmgMesh, mmmgMet );
    }
    else
    {
        if (!muse_metric)
        {
            MMG<TDim>::Init_mesh_ls( mmmgMesh, mmmgLs );
        }
        else
        {
             MMG<TDim>::Init_mesh_ls_met( mmmgMesh, mmmgLs, mmmgMet );
        }
    }
}

template<int TDim>
MMGMesher<TDim>::~MMGMesher()
{
    if (!muse_level_set)
    {
        MMG<TDim>::Free_mesh_met(mmmgMesh, mmmgMet);
    }
    else
    {
        if (!muse_metric)
        {
            MMG<TDim>::Free_mesh_ls(mmmgMesh, mmmgLs);
        }
        else
        {
            MMG<TDim>::Free_mesh_ls_met(mmmgMesh, mmmgLs, mmmgMet);
        }
    }
}

template<int TDim>
void MMGMesher<TDim>::SetValue(const Variable<int>& rVariable, const int& Value)
{
    // TODO
}

template<int TDim>
void MMGMesher<TDim>::SetValue(const int& param, const int& Value)
{
    MMG<TDim>::Set_iparameter(mmmgMesh, mmmgMet, param, Value);
}

template<int TDim>
void MMGMesher<TDim>::SetValue(const Variable<double>& rVariable, const double& Value)
{
    if (rVariable == MMG_GRADATION) // Control the gradation
        MMG<TDim>::Set_dparameter(mmmgMesh, mmmgMet, MMG<TDim>::DPARAM_hgrad, Value);
    else if (rVariable == MMG_HAUSDORFF_DISTANCE) // Control the Hausdorff distance
        MMG<TDim>::Set_dparameter(mmmgMesh, mmmgMet, MMG<TDim>::DPARAM_hausd, Value);
    else if (rVariable == MMG_MINIMAL_MESH_SIZE) // Control the minimal mesh size
        MMG<TDim>::Set_dparameter(mmmgMesh, mmmgMet, MMG<TDim>::DPARAM_hmin, Value);
    else if (rVariable == MMG_MAXIMAL_MESH_SIZE) // Control the minimal mesh size
        MMG<TDim>::Set_dparameter(mmmgMesh, mmmgMet, MMG<TDim>::DPARAM_hmax, Value);
    else if (rVariable == MMG_CONSTANT_MESH_SIZE) // Control the constant mesh size
        MMG<TDim>::Set_dparameter(mmmgMesh, mmmgMet, MMG<TDim>::DPARAM_hsiz, Value);
}

template<int TDim>
void MMGMesher<TDim>::SetValue(const int& param, const double& Value)
{
    MMG<TDim>::Set_dparameter(mmmgMesh, mmmgMet, param, Value);
}

template<int TDim>
void MMGMesher<TDim>::Initialize(ModelPart& r_model_part)
{
    /* Build mesh in MMG5 format */

    /** Manually set of the mesh */
    /** a) give the size of the mesh: vertices, triangles, quads, edges */

    int nvertices = r_model_part.Nodes().size();
    int nelements = r_model_part.Elements().size();
    int cnt;
    int ierr;

    if ( MMG<TDim>::Set_meshSize(mmmgMesh, nvertices, nelements) != 1 )
    {
        KRATOS_ERROR << "MMG failed to initialize by model_part " << r_model_part.Name();
    }

    if (muse_level_set)
    {
        if ( MMG<TDim>::Set_iparameter(mmmgMesh, mmmgLs, MMG<TDim>::IPARAM_iso, 1) != 1 )
            KRATOS_ERROR << "Error setting the level set remeshing mode";

        if ( MMG<TDim>::Set_solSize(mmmgMesh, mmmgLs, MMG5_Vertex, nvertices, MMG5_Scalar) != 1 )
        {
            KRATOS_ERROR << "Error setting the level set vector";
        }
    }

    if (muse_metric)
    {
        /** Manually set of the metric */
        /** give info for the metric structure: sol applied on vertex entities, the sol is scalar*/
        if ( MMG<TDim>::Set_solSize(mmmgMesh, mmmgMet, MMG5_Vertex, nvertices, MMG5_Scalar) != 1 )
        {
            KRATOS_ERROR << "Error setting the refinement metric vector";
        }
    }

    /** give the vertices: for each vertex, give the coordinates, the reference
      and the position in mesh of the vertex */
    cnt = 0;
    for (auto it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
    {
        ++cnt;
        MMG<TDim>::Set_vertex(mmmgMesh, it->X(), it->Y(), it->Z(), 0, cnt);
    }

    std::cout << Info() << ": Read " << cnt << " vertices" << std::endl;

    if (muse_level_set)
    {
        /** give the level set to the vertices */
        cnt = 0;
        for (auto it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
        {
            ++cnt;

            if (it->SolutionStepsDataHas(NODAL_MMG_LEVEL_SET))
                MMG<TDim>::Set_scalarSol(mmmgLs, it->GetSolutionStepValue(NODAL_MMG_LEVEL_SET), cnt);
            else
                KRATOS_ERROR << "The level set is not assigned at node " << it->Id();
        }

        std::cout << Info() << ": Read level set for " << cnt << " vertices" << std::endl;
    }

    if (muse_metric)
    {
        /** give the metric to the vertices */
        cnt = 0;
        for (auto it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
        {
            ++cnt;

            bool has_metric = false;

            if (it->SolutionStepsDataHas(NODAL_MMG_SCALAR_METRIC))
            {
                has_metric = true;
                MMG<TDim>::Set_scalarSol(mmmgMet, it->GetSolutionStepValue(NODAL_MMG_SCALAR_METRIC), cnt);
            }

            if (it->SolutionStepsDataHas(NODAL_MMG_VECTOR_METRIC))
            {
                has_metric = true;
                const Vector& met = it->GetSolutionStepValue(NODAL_MMG_VECTOR_METRIC);
                MMG<TDim>::Set_vectorSol(mmmgMet, met, cnt);
            }

            if (it->SolutionStepsDataHas(NODAL_MMG_TENSOR_METRIC))
            {
                has_metric = true;
                const Matrix& met = it->GetSolutionStepValue(NODAL_MMG_TENSOR_METRIC);
                MMG<TDim>::Set_tensorSol(mmmgMet, met, cnt);
            }

            if (!has_metric)
                KRATOS_ERROR << "The metric is not assigned at node " << it->Id();
        }

        std::cout << Info() << ": Read metric for " << cnt << " vertices" << std::endl;
    }

    /** give the triangle/tetrahedron */
    cnt = 0;
    for (auto it = r_model_part.Elements().begin(); it != r_model_part.Elements().end(); ++it)
    {
        const Properties& rProperties = static_cast<const std::remove_reference_t<decltype(*it)>&>(*it).GetProperties();
        MMG<TDim>::Set_element(mmmgMesh, it->GetGeometry(), static_cast<int>(rProperties.Id()), ++cnt);
    }

    if constexpr (TDim == 2)
        std::cout << Info() << ": Read " << cnt << " triangles" << std::endl;
    else if constexpr (TDim == 3)
        std::cout << Info() << ": Read " << cnt << " tetrahedra" << std::endl;

    cnt = 0;
    for (auto it = r_model_part.Conditions().begin(); it != r_model_part.Conditions().end(); ++it)
    {
        const Properties& rProperties = static_cast<const std::remove_reference_t<decltype(*it)>&>(*it).GetProperties();
        MMG<TDim>::Set_condition(mmmgMesh, it->GetGeometry(), static_cast<int>(rProperties.Id()), ++cnt);
    }

    if constexpr (TDim == 2)
        std::cout << Info() << ": Read " << cnt << " lines" << std::endl;
    else if constexpr (TDim == 3)
        std::cout << Info() << ": Read " << cnt << " triangles" << std::endl;

    if (muse_level_set)
    {
        /** check if the number of given entities match with mesh size */
        if ( MMG<TDim>::Chk_meshData(mmmgMesh, mmmgLs) != 1 )
        {
            KRATOS_ERROR << "Error with MMG mesh";
        }
    }

    if (muse_metric)
    {
        /** check if the number of given entities match with mesh size */
        if ( MMG<TDim>::Chk_meshData(mmmgMesh, mmmgMet) != 1 )
        {
            KRATOS_ERROR << "Error with MMG mesh";
        }
    }

    /** remesh function */
    if (!muse_level_set)
    {
        // WARNING: the MMGxD_mmgxdlib function returns 1 if success, 0 if fail.
        ierr = MMG<TDim>::Remesh(mmmgMesh, mmmgMet);
    }
    else
    {
        ierr = MMG<TDim>::RemeshWithLevelSet(mmmgMesh, mmmgLs, mmmgMet);
    }

    if ( ierr == MMG5_STRONGFAILURE )
    {
        KRATOS_ERROR << "BAD ENDING OF MMGxDLIB: UNABLE TO SAVE MESH";
    }
    else if ( ierr == MMG5_LOWFAILURE )
    {
        KRATOS_ERROR << "BAD ENDING OF MMGxDLIB";
    }

    std::cout << Info() << " is successfully initialized from model_part " << r_model_part.Name() << std::endl;
}

template<int TDim>
void MMGMesher<TDim>::Export(ModelPart& r_model_part, const std::size_t& lastNodeId, const std::size_t& lastElementId,
        const std::string& SampleElementName, const std::size_t& lastPropertiesIndex) const
{
    /// Get the mesh size
    int nvertices, nelements;
    MMG<TDim>::Get_meshSize(mmmgMesh, nvertices, nelements);

    array_1d<double, 3> coords;
    noalias(coords) = ZeroVector(3);
    int ref, isCorner, isRequired;
    std::size_t nodeId = lastNodeId;
    for (int i = 0; i < nvertices; ++i)
    {
        MMG<TDim>::Get_vertex(mmmgMesh, coords, ref, isCorner, isRequired);
        r_model_part.CreateNewNode(++nodeId, coords[0], coords[1], coords[2]);
    }

    std::cout << Info() << ": " << nvertices << " nodes are added to model_part " << r_model_part.Name() << std::endl;

    const Element& SampleElement = KratosComponents<Element>::Get(SampleElementName);

    std::size_t elementId = lastElementId;
    typename Element::NodesArrayType temp_element_nodes;
    std::vector<int> conn(MMG<TDim>::Number_of_nodes_per_element());
    for (int i = 0; i < nelements; ++i)
    {
        MMG<TDim>::Get_element(mmmgMesh, conn, ref, isRequired);

        temp_element_nodes.clear();
        for (int j = 0; j < conn.size(); ++j)
            temp_element_nodes.push_back(r_model_part.pGetNode(conn[j]+lastNodeId));

        Properties::Pointer pProperties = r_model_part.pGetProperties(lastPropertiesIndex+ref);

        Element::Pointer pNewElement = SampleElement.Create(++elementId, temp_element_nodes, pProperties);
        r_model_part.AddElement(pNewElement);
    }

    std::cout << Info() << ": " << nelements << " elements of type " + SampleElementName + " are added to model_part " << r_model_part.Name() << std::endl;
}

template<int TDim>
void MMGMesher<TDim>::SaveMesh(const std::string& Name) const
{
    if (mmmgMesh != NULL)
    {
        if ( MMG<TDim>::saveMesh(mmmgMesh, Name.c_str()) != 1 )
            KRATOS_ERROR << "Error saving mesh to " << Name;
        std::cout << "Saving mesh to " << Name << " successfully" << std::endl;
    }
}

template<int TDim>
void MMGMesher<TDim>::SaveLevelSet(const std::string& Name) const
{
    if (mmmgMesh != NULL && mmmgLs != NULL)
    {
        if ( MMG<TDim>::saveSol(mmmgMesh, mmmgLs, Name.c_str()) != 1 )
            KRATOS_ERROR << "Error saving level set to " << Name;
        std::cout << "Saving level set to " << Name << " successfully" << std::endl;
    }
}

template<int TDim>
void MMGMesher<TDim>::SaveMetric(const std::string& Name) const
{
    if (mmmgMesh != NULL && mmmgMet != NULL)
    {
        if ( MMG<TDim>::saveSol(mmmgMesh, mmmgMet, Name.c_str()) != 1 )
            KRATOS_ERROR << "Error saving metric to " << Name;
        std::cout << "Saving metric to " << Name << " successfully" << std::endl;
    }
}

// template instantiation
template class MMGMesher<2>;
template class MMGMesher<3>;

}
