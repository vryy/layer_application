#ifndef LAYER_APP_MMG_2D_ADAPTER_H_INCLUDED
#define LAYER_APP_MMG_2D_ADAPTER_H_INCLUDED

#include "mmg/mmg2d/libmmg2d.h"

struct MMG2DAdapter
{
    typedef MMG5_pMesh mesh_t;
    typedef MMG5_pSol met_t;
    typedef MMG5_pSol ls_t;
    typedef MMG5_pSol sol_t;

    static constexpr int IPARAM_iso = MMG2D_IPARAM_iso;
    static constexpr int DPARAM_hgrad = MMG2D_DPARAM_hgrad;
    static constexpr int DPARAM_hausd = MMG2D_DPARAM_hausd;
    static constexpr int DPARAM_hmin = MMG2D_DPARAM_hmin;
    static constexpr int DPARAM_hmax = MMG2D_DPARAM_hmax;
    static constexpr int DPARAM_hsiz = MMG2D_DPARAM_hsiz;
    static constexpr int DPARAM_rmc = MMG2D_DPARAM_rmc;

    static inline int Init_mesh_met(mesh_t& p_mesh, met_t& p_met)
    {
        /** ------------------------------------------------------------------ */
        /** Initialisation of mesh and sol structures */
        /* args of InitMesh:
        * MMG5_ARG_start: we start to give the args of a variadic func
        * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
        * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
        * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
        * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */
        return MMG2D_Init_mesh( MMG5_ARG_start,
                         MMG5_ARG_ppMesh, &p_mesh,
                         MMG5_ARG_ppMet, &p_met,
                         MMG5_ARG_end );
    }

    static inline int Free_mesh_met(mesh_t& p_mesh, met_t& p_met)
    {
        /** Free the MMG2D5 structures */
        return MMG2D_Free_all( MMG5_ARG_start,
                        MMG5_ARG_ppMesh, &p_mesh,
                        MMG5_ARG_ppMet, &p_met,
                        MMG5_ARG_end );
    }

    static inline int Free_mesh(mesh_t& p_mesh)
    {
        /** Free the MMG2D5 structures */
        return MMG2D_Free_all( MMG5_ARG_start,
                        MMG5_ARG_ppMesh, &p_mesh,
                        MMG5_ARG_end );
    }

    static inline int Init_mesh_ls(mesh_t& p_mesh, ls_t& p_ls)
    {
        /** ------------------------------------------------------------------ */
        /** Initialisation of mesh and sol structures */
        /* args of InitMesh:
        * MMG5_ARG_start: we start to give the args of a variadic func
        * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
        * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
        * MMG5_ARG_ppLs: next arg will be a pointer over a MMG5_pSol storing a level set
        * &mmgLs: pointer toward your MMG5_pSol (that store your level set) */
        return MMG2D_Init_mesh( MMG5_ARG_start,
                         MMG5_ARG_ppMesh, &p_mesh,
                         MMG5_ARG_ppLs, &p_ls,
                         MMG5_ARG_end );
    }

    static inline int Free_mesh_ls(mesh_t& p_mesh, ls_t& p_ls)
    {
        /** Free the MMG2D5 structures */
        return MMG2D_Free_all( MMG5_ARG_start,
                        MMG5_ARG_ppMesh, &p_mesh,
                        MMG5_ARG_ppLs, &p_ls,
                        MMG5_ARG_end );
    }

    static inline int Init_mesh_ls_met(mesh_t& p_mesh, ls_t& p_ls, met_t& p_met)
    {
        /** ------------------------------------------------------------------ */
        /** Initialisation of mesh and sol structures */
        /* args of InitMesh:
        * MMG5_ARG_start: we start to give the args of a variadic func
        * MMG5_ARG_ppMesh: next arg will be a pointer over a MMG5_pMesh
        * &mmgMesh: pointer toward your MMG5_pMesh (that store your mesh)
        * MMG5_ARG_ppLs: next arg will be a pointer over a MMG5_pSol storing a level set
        * &mmgLs: pointer toward your MMG5_pSol (that store your level set)
        * MMG5_ARG_ppMet: next arg will be a pointer over a MMG5_pSol storing a metric
        * &mmgSol: pointer toward your MMG5_pSol (that store your metric) */
        return MMG2D_Init_mesh( MMG5_ARG_start,
                         MMG5_ARG_ppMesh, &p_mesh,
                         MMG5_ARG_ppLs, &p_ls,
                         MMG5_ARG_ppMet, &p_met,
                         MMG5_ARG_end );
    }

    static inline int Free_mesh_ls_met(mesh_t& p_mesh, ls_t& p_ls, met_t& p_met)
    {
        /** Free the MMG2D5 structures */
        return MMG2D_Free_all( MMG5_ARG_start,
                        MMG5_ARG_ppMesh, &p_mesh,
                        MMG5_ARG_ppLs, &p_ls,
                        MMG5_ARG_ppMet, &p_met,
                        MMG5_ARG_end );
    }

    static inline int Set_iparameter(mesh_t& p_mesh, sol_t& p_sol, const int& iparam, const int& val)
    {
        return MMG2D_Set_iparameter(p_mesh, p_sol, iparam, val);
    }

    static inline int Set_dparameter(mesh_t& p_mesh, sol_t& p_sol, const int& dparam, const double& val)
    {
        return MMG2D_Set_dparameter(p_mesh, p_sol, dparam, val);
    }

    static inline int Set_meshSize(mesh_t& p_mesh, const int& nvertices, const int& nelements)
    {
        const int nquads = 0;
        const int nedges = 0;
        return MMG2D_Set_meshSize(p_mesh, nvertices, nelements, nquads, nedges);
    }

    static inline int Set_meshSize(mesh_t& p_mesh, const int& nvertices, const int& nelements, const int& nconditions)
    {
        const int nquads = 0;
        return MMG2D_Set_meshSize(p_mesh, nvertices, nelements, nquads, nconditions);
    }

    static inline int Get_meshSize(const mesh_t& p_mesh, int& nvertices, int& nelements)
    {
        int nquads, nedges;
        return MMG2D_Get_meshSize(p_mesh, &nvertices, &nelements, &nquads, &nedges);
    }

    static inline int Get_meshSize(const mesh_t& p_mesh, int& nvertices, int& nelements, int& nconditions)
    {
        int nquads;
        return MMG2D_Get_meshSize(p_mesh, &nvertices, &nelements, &nquads, &nconditions);
    }

    static inline int Set_solSize(mesh_t& p_mesh, sol_t& p_sol, const int& typEntity, const int& nvertices, const int& typSol)
    {
        return MMG2D_Set_solSize(p_mesh, p_sol, typEntity, nvertices, typSol);
    }

    static inline int Set_vertex(mesh_t& p_mesh, const double& x, const double& y, const double& z, const int& ref, const int& pos)
    {
        return MMG2D_Set_vertex(p_mesh, x, y, ref, pos);
    }

    template<typename TVectorType>
    static inline int Get_vertex(const mesh_t& p_mesh, TVectorType& coords, int& ref, int& isCorner, int& isRequired)
    {
        return MMG2D_Get_vertex(p_mesh, &coords[0], &coords[1], &ref, &isCorner, &isRequired);
    }

    static inline int Set_scalarSol(sol_t& p_sol, const double& s, const int& pos)
    {
        return MMG2D_Set_scalarSol(p_sol, s, pos);
    }

    template<typename TVectorType>
    static inline int Set_vectorSol(sol_t& p_sol, const TVectorType& vec, const int& pos)
    {
        return MMG2D_Set_vectorSol(p_sol, vec[0], vec[1], pos);
    }

    template<typename TTensorType>
    static inline int Set_tensorSol(sol_t& p_sol, const TTensorType& ten, const int& pos)
    {
        return MMG2D_Set_tensorSol(p_sol, ten(0, 0), ten(0, 1), ten(1, 1), pos);
    }

    template<typename TGeometryType>
    static inline int Set_element(mesh_t& p_mesh, const TGeometryType& r_geom, const int& ref, const int& pos)
    {
        return MMG2D_Set_triangle(p_mesh, static_cast<int>(r_geom[0].Id()), static_cast<int>(r_geom[1].Id()), static_cast<int>(r_geom[2].Id()), ref, pos);
    }

    template<typename TGeometryType>
    static inline int Set_condition(mesh_t& p_mesh, const TGeometryType& r_geom, const int& ref, const int& pos)
    {
        return MMG2D_Set_edge(p_mesh, static_cast<int>(r_geom[0].Id()), static_cast<int>(r_geom[1].Id()), ref, pos);
    }

    static inline int Set_multiMat(mesh_t& p_mesh, sol_t& p_sol, const int& ref, const int& split, const int& rin, const int& rex)
    {
        return MMG2D_Set_multiMat(p_mesh, p_sol, ref, split, rin, rex);
    }

    template<typename TConnectivityType>
    static inline int Get_element(const mesh_t& p_mesh, TConnectivityType& conn, int& ref, int& isRequired)
    {
        return MMG2D_Get_triangle(p_mesh, &conn[0], &conn[1], &conn[2], &ref, &isRequired);
    }

    static inline int Get_elements(const mesh_t& p_mesh, int* conns, int* refs, int* areRequired)
    {
        return MMG2D_Get_triangles(p_mesh, conns, refs, areRequired);
    }

    template<typename TConnectivityType>
    static inline int Get_condition(const mesh_t& p_mesh, TConnectivityType& conn, int& ref, int& isRequired)
    {
        int isRidge;
        return MMG2D_Get_edge(p_mesh, &conn[0], &conn[1], &ref, &isRidge, &isRequired);
    }

    static inline int Get_parent_element(mesh_t& p_mesh, const int& cond_id, int& parent_elem_id, int& local_index)
    {
        return MMG2D_Get_triFromEdge(p_mesh, cond_id, &parent_elem_id, &local_index);
    }

    static inline int Get_conditions(const mesh_t& p_mesh, int* conns, int* refs, int* areRequired)
    {
        int* areRidges;
        return MMG2D_Get_edges(p_mesh, conns, refs, areRidges, areRequired);
    }

    static inline constexpr int Number_of_nodes_per_element()
    {
        return 3;
    }

    static inline constexpr int Number_of_nodes_per_condition()
    {
        return 2;
    }

    static inline int Chk_meshData(mesh_t& p_mesh, sol_t& p_sol)
    {
        return MMG2D_Chk_meshData(p_mesh, p_sol);
    }

    static inline int Remesh(mesh_t& p_mesh, met_t& p_met)
    {
        return MMG2D_mmg2dlib(p_mesh, p_met);
    }

    static inline int RemeshWithLevelSet(mesh_t& p_mesh, ls_t& p_ls, met_t& p_met)
    {
        return MMG2D_mmg2dls(p_mesh, p_ls, p_met);
    }

    static inline int saveMesh(const mesh_t& p_mesh, const char* filename)
    {
        return MMG2D_saveMesh(p_mesh, filename);
    }

    static inline int saveSol(const mesh_t& p_mesh, const sol_t& p_sol, const char *filename)
    {
        return MMG2D_saveSol(p_mesh, p_sol, filename);
    }

    static inline int saveMsh(const mesh_t& p_mesh, const sol_t& p_sol, const char* filename)
    {
        return MMG2D_saveMshMesh(p_mesh, p_sol, filename);
    }
};

#endif // LAYER_APP_MMG_2D_ADAPTER_H_INCLUDED
