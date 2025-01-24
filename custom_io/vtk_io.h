//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         layer_application/LICENSE.txt
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//


#if !defined(KRATOS_VTK_IO_H_INCLUDED)
#define  KRATOS_VTK_IO_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <iomanip>

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"
#include "utilities/timer.h"
#include "custom_utilities/vtk.h"

namespace Kratos
{

/**
 * This class defines an interface to Paraview
 * in order to provide VTK compliant I/O functionality
 * 14 Aug 2017: supports only nodal scalar & vector results
 */
template<class TMeshContainer>
class VtkIO // : public IO
{
public:
    ///pointer definition of VtkIO
    KRATOS_CLASS_POINTER_DEFINITION(VtkIO);

    ///typedefs
//    typedef IO BaseType;
    typedef ModelPart::MeshType MeshType;
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    typedef ModelPart::NodesContainerType NodesContainerType;
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;
    typedef GeometryData::IntegrationMethod IntegrationMethodType;
    typedef GeometryData::KratosGeometryFamily KratosGeometryFamily;

    /// Constructor
    VtkIO( const std::string& rDatafilename, const VTK_PostMode& rMode )
    {
        mMode = rMode;
        mResultFileOpen = false;
        mResultFileName = rDatafilename;
        mResultCellProperties = false;
        InitializeResultFile(mResultFileName);
        SetUpMeshContainers();
    }

    /// Destructor.
    virtual ~VtkIO()
    {
        Timer::PrintTimingInformation();

        this->CloseResultFile();
    }

    ///initialization functions
    /**
     * creates the mesh containers for all different element types.
     * Note that the containers are not filled yet in here!
     */
    void SetUpMeshContainers()
    {
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,
                                          VTK_Hexahedron, "Kratos_Hexahedra3D20_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D27,
                                          VTK_Hexahedron, "Kratos_Hexahedra3D27_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,
                                          VTK_Hexahedron, "Kratos_Hexahedra3D8_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Prism3D15,
                                          VTK_Wedge, "Kratos_Prism3D15_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Prism3D6,
                                          VTK_Wedge, "Kratos_Prism3D6_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4,
                                          VTK_Quad, "Kratos_Quadrilateral2D4_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8,
                                          VTK_Quad, "Kratos_Quadrilateral2D8_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9,
                                          VTK_Quad, "Kratos_Quadrilateral2D9_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4,
                                          VTK_Quad, "Kratos_Quadrilateral3D4_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8,
                                          VTK_Quad, "Kratos_Quadrilateral3D8_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9,
                                          VTK_Quad, "Kratos_Quadrilateral3D9_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,
                                          VTK_Tetra, "Kratos_Tetrahedra3D10_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,
                                          VTK_Tetra, "Kratos_Tetrahedra3D4_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle2D3,
                                          VTK_Triangle, "Kratos_Triangle2D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle2D6,
                                          VTK_Triangle, "Kratos_Triangle2D6_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle3D3,
                                          VTK_Triangle, "Kratos_Triangle3D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triangle3D6,
                                          VTK_Triangle, "Kratos_Triangle3D6_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line2D2,
                                          VTK_Line, "Kratos_Line2D2_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line3D2,
                                          VTK_Line, "Kratos_Line3D2_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line2D3,
                                          VTK_Line, "Kratos_Line2D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Line3D3,
                                          VTK_Line, "Kratos_Line3D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Point3D,
                                          VTK_Pixel, "Kratos_Point3D_Mesh" ) );
        #ifndef SD_APP_FORWARD_COMPATIBILITY
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tritagon,
                                          VTK_Polygon, "Kratos_Tritagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tetragon,
                                          VTK_Polygon, "Kratos_Tetragon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Pentagon,
                                          VTK_Polygon, "Kratos_Pentagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexagon,
                                          VTK_Polygon, "Kratos_Hexagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Heptagon,
                                          VTK_Polygon, "Kratos_Heptagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Octagon,
                                          VTK_Polygon, "Kratos_Octagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Nonagon,
                                          VTK_Polygon, "Kratos_Nonagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Decagon,
                                          VTK_Polygon, "Kratos_Decagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hendecagon,
                                          VTK_Polygon, "Kratos_Hendecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Dodecagon,
                                          VTK_Polygon, "Kratos_Dodecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Triskaidecagon,
                                          VTK_Polygon, "Kratos_Triskaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Tetrakaidecagon,
                                          VTK_Polygon, "Kratos_Tetrakaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Pentadecagon,
                                          VTK_Polygon, "Kratos_Pentadecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Hexakaidecagon,
                                          VTK_Polygon, "Kratos_Hexakaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Heptadecagon,
                                          VTK_Polygon, "Kratos_Heptadecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Octakaidecagon,
                                          VTK_Polygon, "Kratos_Octakaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Enneadecagon,
                                          VTK_Polygon, "Kratos_Enneadecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Octahedron,
                                          VTK_Polyhedron, "Kratos_Octahedron_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Dodecahedron,
                                          VTK_Polyhedron, "Kratos_Dodecahedron_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::KratosGeometryType::Kratos_Icosahedron,
                                          VTK_Polyhedron, "Kratos_Icosahedron_Mesh" ) );
        #endif
    } // SetUpMeshContainers

    /**
     * sets up the file name
     */
    void InitializeResultFile( std::string const& rResultFileName )
    {
        //std::cout << "initializing result files" << std::endl;
        mResultFileName = rResultFileName;
    }

    /**
     * close the file
     */
    void CloseResultFile()
    {
        if ( mResultFileOpen )
        {
            fclose(mResultFile);
            mResultFileOpen = false;
        }
    }

    /**
     * Turn back information as a string.
     */
    virtual std::string Info() const
    {
        return "VTK IO";
    }

    /**
     * Print information about this object.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Number of vtk mesh container: " << mVtkMeshContainers.size();
    }

    /**
     * This has to be called for each solution step BEFORE any results
     * (on nodes) is written
     * @param name the current solution step (i.e. time)
     * @param rThisMesh the mesh containing the results
     */
    virtual void Initialize( double name, const MeshType& rThisMesh )
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /**
     * This has to be called for each solution step to write all the results to file
     * have been written
     */
    virtual void Finalize()
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /**
     * Register the nodal results (to be written later in Finalize)
     */
    void RegisterNodalResults(const Variable<double>& rVariable, const std::size_t SolutionStepNumber)
    {
        for (std::size_t i = 0; i < mNodalDoubleVariablesMap.size(); ++i)
        {
            if (mNodalDoubleVariablesMap[i].first == rVariable
             && mNodalDoubleVariablesMap[i].second == SolutionStepNumber)
            {
                std::cout << "Info: Nodal result " << rVariable.Name() << " at solution step " << SolutionStepNumber
                          << " is registerred" << std::endl;
                return;
            }
        }

        mNodalDoubleVariablesMap.push_back(std::make_pair(rVariable, SolutionStepNumber));
    }

    /**
     * Register the nodal results (to be written later in Finalize)
     */
    void RegisterNodalResults(const Variable<array_1d<double, 3> >& rVariable, const std::size_t SolutionStepNumber)
    {
        for (std::size_t i = 0; i < mNodalArray1dVariablesMap.size(); ++i)
        {
            if (mNodalArray1dVariablesMap[i].first == rVariable
             && mNodalArray1dVariablesMap[i].second == SolutionStepNumber)
            {
                std::cout << "Info: Nodal result " << rVariable.Name() << " at solution step " << SolutionStepNumber
                          << " is registerred" << std::endl;
                return;
            }
        }

        mNodalArray1dVariablesMap.push_back(std::make_pair(rVariable, SolutionStepNumber));
    }

    /**
     * Register the nodal results (to be written later in Finalize)
     */
    void RegisterNodalResults(const Variable<Vector>& rVariable, const std::size_t SolutionStepNumber, const std::size_t vec_size)
    {
        for (std::size_t i = 0; i < mNodalVectorVariablesMap.size(); ++i)
        {
            if (mNodalVectorVariablesMap[i].first == rVariable
             && mNodalVectorVariablesMap[i].second == SolutionStepNumber)
            {
                std::cout << "Info: Nodal result " << rVariable.Name() << " at solution step " << SolutionStepNumber
                          << " is registerred" << std::endl;
                return;
            }
        }

        mNodalVectorVariablesMap.push_back(std::make_pair(rVariable, SolutionStepNumber));
        mNodalVectorVariablesSizeMap.push_back(std::make_pair(rVariable, vec_size));
    }

    /**
     * Register the cell results (to be written later in Finalize)
     */
    void RegisterCellResults(const Variable<int>& rVariable, const std::size_t SolutionStepNumber)
    {
        for (std::size_t i = 0; i < mCellIntegerVariablesMap.size(); ++i)
        {
            if (mCellIntegerVariablesMap[i].first == rVariable
             && mCellIntegerVariablesMap[i].second == SolutionStepNumber)
            {
                std::cout << "Info: Cell result " << rVariable.Name() << " at solution step " << SolutionStepNumber
                          << " is registerred" << std::endl;
                return;
            }
        }

        mCellIntegerVariablesMap.push_back(std::make_pair(rVariable, SolutionStepNumber));
    }

    /**
     * Register the cell results (to be written later in Finalize)
     */
    void RegisterCellResults(const Variable<double>& rVariable, const std::size_t SolutionStepNumber)
    {
        for (std::size_t i = 0; i < mCellDoubleVariablesMap.size(); ++i)
        {
            if (mCellDoubleVariablesMap[i].first == rVariable
             && mCellDoubleVariablesMap[i].second == SolutionStepNumber)
            {
                std::cout << "Info: Cell result " << rVariable.Name() << " at solution step " << SolutionStepNumber
                          << " is registerred" << std::endl;
                return;
            }
        }

        mCellDoubleVariablesMap.push_back(std::make_pair(rVariable, SolutionStepNumber));
    }

    /**
     * Register to write the cell Properties (to be written later in Finalize)
     */
    void RegisterCellProperties()
    {
        mResultCellProperties = true;
    }

protected:

    bool mResultFileOpen;
    std::string mResultFileName; // result file name
    FILE* mResultFile; // pointer to file writer
    VTK_PostMode mMode; // write flag

    std::vector<TMeshContainer>& MeshContainers() {return mVtkMeshContainers;}

    void BeginNodalResultsHeader( FILE* pResultFile ) const
    {
        std::stringstream ScalarVars;
        std::stringstream VectorVars;

        for(auto it = mNodalDoubleVariablesMap.begin(); it != mNodalDoubleVariablesMap.end(); ++it)
        {
            ScalarVars << it->first.Name() << ",";
        }

        for(auto it = mNodalArray1dVariablesMap.begin(); it != mNodalArray1dVariablesMap.end(); ++it)
        {
            VectorVars << it->first.Name() << ",";
        }

        for(auto it = mNodalVectorVariablesMap.begin(); it != mNodalVectorVariablesMap.end(); ++it)
        {
            VectorVars << it->first.Name() << ",";
        }

        fprintf( pResultFile, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n", ScalarVars.str().c_str(), VectorVars.str().c_str() );
    }

    void EndNodalResultsHeader( FILE* pResultFile ) const
    {
        fprintf( pResultFile, "      </PointData>\n" );
    }

    void BeginCellResultsHeader( FILE* pResultFile ) const
    {
        std::stringstream ScalarVars;
        std::stringstream VectorVars;

        for(auto it = mCellIntegerVariablesMap.begin(); it != mCellIntegerVariablesMap.end(); ++it)
        {
            ScalarVars << it->first.Name() << ",";
        }

        for(auto it = mCellDoubleVariablesMap.begin(); it != mCellDoubleVariablesMap.end(); ++it)
        {
            ScalarVars << it->first.Name() << ",";
        }

        fprintf( pResultFile, "      <CellData Scalars=\"%s\" Vectors=\"%s\">\n", ScalarVars.str().c_str(), VectorVars.str().c_str() );
    }

    void EndCellResultsHeader( FILE* pResultFile ) const
    {
        fprintf( pResultFile, "      </CellData>\n" );
    }

    /**
     * writes nodal results for all variables at once
     */
    void WriteNodalResults( FILE* pResultFile, const NodesContainerType& rNodes ) const
    {
        for(auto it = mNodalDoubleVariablesMap.begin(); it != mNodalDoubleVariablesMap.end(); ++it)
        {
            WriteNodalResult( pResultFile, it->first, rNodes, it->second );
        }

        for(auto it = mNodalArray1dVariablesMap.begin(); it != mNodalArray1dVariablesMap.end(); ++it)
        {
            WriteNodalResult( pResultFile, it->first, rNodes, it->second );
        }

        std::size_t cnt = 0;
        for(auto it = mNodalVectorVariablesMap.begin(); it != mNodalVectorVariablesMap.end(); ++it, ++cnt)
        {
            std::size_t size = mNodalVectorVariablesSizeMap[cnt].second;
            WriteNodalResult( pResultFile, it->first, rNodes, it->second, size );
        }
    }

    /**
     * writes nodal results for all variables at once
     */
    template<typename EntitiesContainerType>
    void WriteCellResults( FILE* pResultFile, const EntitiesContainerType& rElements ) const
    {
        for(auto it = mCellIntegerVariablesMap.begin(); it != mCellIntegerVariablesMap.end(); ++it)
        {
            WriteCellResult<int>( pResultFile, it->first, rElements, it->second );
        }

        for(auto it = mCellDoubleVariablesMap.begin(); it != mCellDoubleVariablesMap.end(); ++it)
        {
            WriteCellResult<double>( pResultFile, it->first, rElements, it->second );
        }

        if (mResultCellProperties)
            WriteCellProperties( pResultFile, rElements, 0 );
    }

private:

    /**
     * member variables
     */
    std::vector<TMeshContainer> mVtkMeshContainers; // each mesh container will go into a piece

    std::vector< std::pair<Variable<double>, std::size_t> > mNodalDoubleVariablesMap;
    std::vector< std::pair<Variable<array_1d<double, 3> >, std::size_t> > mNodalArray1dVariablesMap;
    std::vector< std::pair<Variable<Vector>, std::size_t> > mNodalVectorVariablesMap;
    std::vector< std::pair<Variable<Vector>, std::size_t> > mNodalVectorVariablesSizeMap;

    std::vector< std::pair<Variable<int>, std::size_t> > mCellIntegerVariablesMap;
    std::vector< std::pair<Variable<double>, std::size_t> > mCellDoubleVariablesMap;

    bool mResultCellProperties; // special flag to write the cell properties

    /**
     * assignment operator
     */
    VtkIO& operator=(VtkIO const& rOther);

    /**
     * writes single nodal result for variable of type double
     */
    void WriteNodalResult( FILE* pResultFile,
                           Variable<double> const& rVariable,
                           const NodesContainerType& rNodes,
                           const std::size_t SolutionStepNumber ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), 1, mMode );

        if (mMode == VTK_PostAscii)
        {
            for ( NodesContainerType::const_iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node )
                VTK_fWriteScalar( pResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable, SolutionStepNumber) );
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(rNodes.end() - rNodes.begin());
            for ( NodesContainerType::const_iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node )
                data_list.push_back(i_node->GetSolutionStepValue(rVariable, SolutionStepNumber));
            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( pResultFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }

        VTK_fEndResult( pResultFile );

        Timer::Stop("Writing Results");
    }

    /**
     * writes single nodal result for variable of type array_1d<double, 3>
     */
    void WriteNodalResult( FILE* pResultFile,
                           Variable<array_1d<double, 3> > const& rVariable,
                           const NodesContainerType& rNodes,
                           const std::size_t SolutionStepNumber ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), 3, mMode );

        if (mMode == VTK_PostAscii)
        {
            for ( NodesContainerType::const_iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node )
            {
                array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable, SolutionStepNumber );
                VTK_fWriteVector3( pResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
            }
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(3*(rNodes.end() - rNodes.begin()));
            for ( NodesContainerType::const_iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node )
            {
                array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable, SolutionStepNumber );
                data_list.push_back(temp[0]);
                data_list.push_back(temp[1]);
                data_list.push_back(temp[2]);
            }
            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( pResultFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }

        VTK_fEndResult( pResultFile );

        Timer::Stop("Writing Results");
    }

    /**
     * writes single nodal result for variable of type Vector
     */
    void WriteNodalResult( FILE* pResultFile,
                           Variable<Vector> const& rVariable,
                           const NodesContainerType& rNodes,
                           const std::size_t SolutionStepNumber,
                           const std::size_t NumberOfComponents ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), NumberOfComponents, mMode );

        if (mMode == VTK_PostAscii)
        {
            double* temp = new double[NumberOfComponents];
            for ( NodesContainerType::const_iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node )
            {
                const Vector& solution = i_node->GetSolutionStepValue( rVariable, SolutionStepNumber );
                for ( unsigned int i = 0; i < NumberOfComponents; ++i )
                    temp[i] = solution(i);
                VTK_fWriteVector( pResultFile, i_node->Id(), NumberOfComponents, temp );
            }
            delete temp;
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(NumberOfComponents*(rNodes.end() - rNodes.begin()));
            for ( NodesContainerType::const_iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node )
            {
                const Vector& solution = i_node->GetSolutionStepValue( rVariable, SolutionStepNumber );
                for ( unsigned int i = 0; i < NumberOfComponents; ++i )
                    data_list.push_back(solution(i));
            }
            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( pResultFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }

        VTK_fEndResult( pResultFile );

        Timer::Stop("Writing Results");
    }

    /**
     * writes single cell result for variable of type scalar (bool|int|double)
     */
    template<typename TScalarType, typename EntitiesContainerType>
    void WriteCellResult( FILE* pResultFile,
                           Variable<TScalarType> const& rVariable,
                           const EntitiesContainerType& rElements,
                           const std::size_t SolutionStepNumber ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), 1, mMode );

        if (mMode == VTK_PostAscii)
        {
            for ( typename EntitiesContainerType::const_iterator i_elem = rElements.begin(); i_elem != rElements.end() ; ++i_elem )
                VTK_fWriteScalar( pResultFile, i_elem->Id(), i_elem->GetValue(rVariable) );
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(rElements.end() - rElements.begin());
            for ( typename EntitiesContainerType::const_iterator i_elem = rElements.begin(); i_elem != rElements.end() ; ++i_elem )
                data_list.push_back(i_elem->GetValue(rVariable));
            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( pResultFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }

        VTK_fEndResult( pResultFile );

        Timer::Stop("Writing Results");
    }

    /**
     * writes cell Properties
     */
    template<typename EntitiesContainerType>
    void WriteCellProperties( FILE* pResultFile,
                              const EntitiesContainerType& rElements,
                              const std::size_t SolutionStepNumber ) const
    {
        Timer::Start("Writing Results");

        const std::string var_name = "Properties";
        VTK_fBeginResult( pResultFile, (char*)(var_name.c_str()), 1, mMode );

        if (mMode == VTK_PostAscii)
        {
            for ( typename EntitiesContainerType::const_iterator i_elem = rElements.begin(); i_elem != rElements.end() ; ++i_elem )
                VTK_fWriteScalar( pResultFile, i_elem->Id(), i_elem->GetProperties().Id() );
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(rElements.end() - rElements.begin());
            for ( typename EntitiesContainerType::const_iterator i_elem = rElements.begin(); i_elem != rElements.end() ; ++i_elem )
                data_list.push_back(i_elem->GetProperties().Id());
            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( pResultFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }

        VTK_fEndResult( pResultFile );

        Timer::Stop("Writing Results");
    }

    /**
     * Copy constructor
     */
    VtkIO(VtkIO const& rOther);

}; // Class VtkIO

/**
 * Input and output
 */

/**
 * output stream function
 */
template<class TMeshContainer>
inline std::ostream& operator << (std::ostream& rOStream, const VtkIO<TMeshContainer>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_VTK_IO_H_INCLUDED  defined
