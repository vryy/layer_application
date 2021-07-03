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
#include "includes/io.h"
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
class VtkIO// : public IO
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
    VtkIO( const std::string& rDatafilename, const VTK_PostMode& rMode)
    {
        mMode = rMode;
        mResultFileOpen = false;
        mResultFileName = rDatafilename;
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
                                          GeometryData::Kratos_Hexahedra3D20,
                                          VTK_Hexahedron, "Kratos_Hexahedra3D20_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexahedra3D27,
                                          VTK_Hexahedron, "Kratos_Hexahedra3D27_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexahedra3D8,
                                          VTK_Hexahedron, "Kratos_Hexahedra3D8_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Prism3D15,
                                          VTK_Wedge, "Kratos_Prism3D15_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Prism3D6,
                                          VTK_Wedge, "Kratos_Prism3D6_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral2D4,
                                          VTK_Quad, "Kratos_Quadrilateral2D4_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral2D8,
                                          VTK_Quad, "Kratos_Quadrilateral2D8_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral2D9,
                                          VTK_Quad, "Kratos_Quadrilateral2D9_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral3D4,
                                          VTK_Quad, "Kratos_Quadrilateral3D4_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral3D8,
                                          VTK_Quad, "Kratos_Quadrilateral3D8_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Quadrilateral3D9,
                                          VTK_Quad, "Kratos_Quadrilateral3D9_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tetrahedra3D10,
                                          VTK_Tetra, "Kratos_Tetrahedra3D10_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tetrahedra3D4,
                                          VTK_Tetra, "Kratos_Tetrahedra3D4_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle2D3,
                                          VTK_Triangle, "Kratos_Triangle2D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle2D6,
                                          VTK_Triangle, "Kratos_Triangle2D6_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle3D3,
                                          VTK_Triangle, "Kratos_Triangle3D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triangle3D6,
                                          VTK_Triangle, "Kratos_Triangle3D6_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line2D2,
                                          VTK_Line, "Kratos_Line2D2_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line3D2,
                                          VTK_Line, "Kratos_Line3D2_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line2D3,
                                          VTK_Line, "Kratos_Line2D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Line3D3,
                                          VTK_Line, "Kratos_Line3D3_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Point3D,
                                          VTK_Pixel, "Kratos_Point3D_Mesh" ) );
        #ifndef SD_APP_FORWARD_COMPATIBILITY
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tritagon,
                                          VTK_Polygon, "Kratos_Tritagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tetragon,
                                          VTK_Polygon, "Kratos_Tetragon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Pentagon,
                                          VTK_Polygon, "Kratos_Pentagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexagon,
                                          VTK_Polygon, "Kratos_Hexagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Heptagon,
                                          VTK_Polygon, "Kratos_Heptagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Octagon,
                                          VTK_Polygon, "Kratos_Octagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Nonagon,
                                          VTK_Polygon, "Kratos_Nonagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Decagon,
                                          VTK_Polygon, "Kratos_Decagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hendecagon,
                                          VTK_Polygon, "Kratos_Hendecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Dodecagon,
                                          VTK_Polygon, "Kratos_Dodecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Triskaidecagon,
                                          VTK_Polygon, "Kratos_Triskaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Tetrakaidecagon,
                                          VTK_Polygon, "Kratos_Tetrakaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Pentadecagon,
                                          VTK_Polygon, "Kratos_Pentadecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Hexakaidecagon,
                                          VTK_Polygon, "Kratos_Hexakaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Heptadecagon,
                                          VTK_Polygon, "Kratos_Heptadecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Octakaidecagon,
                                          VTK_Polygon, "Kratos_Octakaidecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Enneadecagon,
                                          VTK_Polygon, "Kratos_Enneadecagon_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Octahedron,
                                          VTK_Polyhedron, "Kratos_Octahedron_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Dodecahedron,
                                          VTK_Polyhedron, "Kratos_Dodecahedron_Mesh" ) );
        mVtkMeshContainers.push_back( TMeshContainer(
                                          GeometryData::Kratos_Icosahedron,
                                          VTK_Polyhedron, "Kratos_Icosahedron_Mesh" ) );
        #endif
    }//SetUpMeshContainers

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
    virtual void Initialize( double name, MeshType& rThisMesh )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /**
     * This has to be called for each solution step to write all the results to file
     * have been written
     */
    virtual void Finalize()
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class function", __FUNCTION__)
    }

    /**
     * Register the nodal results (to be written later in Finalize)
     */
    void RegisterNodalResults(Variable<double> const& rVariable, const std::size_t& SolutionStepNumber)
    {
        mNodalDoubleVariablesMap[rVariable] = SolutionStepNumber;
    }

    /**
     * Register the nodal results (to be written later in Finalize)
     */
    void RegisterNodalResults(Variable<array_1d<double, 3> > const& rVariable, const std::size_t& SolutionStepNumber)
    {
        mNodalArray1dVariablesMap[rVariable] = SolutionStepNumber;
    }

    /**
     * Register the nodal results (to be written later in Finalize)
     */
    void RegisterNodalResults(Variable<Vector> const& rVariable, const std::size_t& SolutionStepNumber, const std::size_t& vec_size)
    {
        mNodalVectorVariablesMap[rVariable] = SolutionStepNumber;
        mNodalVectorVariablesSizeMap[rVariable] = vec_size;
    }

protected:
    /**
     * File names
     */
    std::string mResultFileName;

    FILE* mResultFile;

    /**
     * Flags
     */
    VTK_PostMode mMode;

    /**
     * member variables
     */
    std::vector<TMeshContainer> mVtkMeshContainers; // each mesh container will go into a piece
    bool mResultFileOpen;

    std::map<Variable<double>, std::size_t> mNodalDoubleVariablesMap;
    std::map<Variable<array_1d<double, 3> >, std::size_t> mNodalArray1dVariablesMap;
    std::map<Variable<Vector>, std::size_t> mNodalVectorVariablesMap;
    std::map<Variable<Vector>, std::size_t> mNodalVectorVariablesSizeMap;

    void BeginResultsHeader( FILE* pResultFile ) const
    {
        std::stringstream ScalarVars;
        std::stringstream VectorVars;

        for(std::map<Variable<double>, std::size_t>::const_iterator it = mNodalDoubleVariablesMap.begin();
                it != mNodalDoubleVariablesMap.end(); ++it)
        {
            ScalarVars << it->first.Name() << ",";
        }

        for(std::map<Variable<array_1d<double, 3> >, std::size_t>::const_iterator it = mNodalArray1dVariablesMap.begin();
                it != mNodalArray1dVariablesMap.end(); ++it)
        {
            VectorVars << it->first.Name() << ",";
        }

        for(std::map<Variable<Vector>, std::size_t>::const_iterator it = mNodalVectorVariablesMap.begin();
                it != mNodalVectorVariablesMap.end(); ++it)
        {
            VectorVars << it->first.Name() << ",";
        }

        fprintf( pResultFile, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n", ScalarVars.str().c_str(), VectorVars.str().c_str() );
    }

    void EndResultsHeader( FILE* pResultFile ) const
    {
        fprintf( pResultFile, "      </PointData>\n" );
    }

    /**
     * writes nodal results for all variables at once
     */
    void WriteNodalResults( FILE* pResultFile, NodesContainerType& rNodes ) const
    {
        for(std::map<Variable<double>, std::size_t>::const_iterator it = mNodalDoubleVariablesMap.begin();
                it != mNodalDoubleVariablesMap.end(); ++it)
        {
            WriteNodalResults( pResultFile, it->first, rNodes, it->second );
        }

        for(std::map<Variable<array_1d<double, 3> >, std::size_t>::const_iterator it = mNodalArray1dVariablesMap.begin();
                it != mNodalArray1dVariablesMap.end(); ++it)
        {
            WriteNodalResults( pResultFile, it->first, rNodes, it->second );
        }

        for(std::map<Variable<Vector>, std::size_t>::const_iterator it = mNodalVectorVariablesMap.begin();
                it != mNodalVectorVariablesMap.end(); ++it)
        {
            auto it_size = mNodalVectorVariablesSizeMap.find(it->first);
            WriteNodalResults( pResultFile, it->first, rNodes, it->second, it_size->second );
        }
    }

    /**
     * writes nodal results for variables of type double
     */
    void WriteNodalResults( FILE* pResultFile,
                            Variable<double> const& rVariable,
                            NodesContainerType& rNodes,
                            const std::size_t& SolutionStepNumber ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), 1, mMode );

        if (mMode == VTK_PostAscii)
        {
            for ( NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node)
                VTK_fWriteScalar( pResultFile, i_node->Id(), i_node->GetSolutionStepValue(rVariable, SolutionStepNumber) );
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(rNodes.end() - rNodes.begin());
            for ( NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node)
                data_list.push_back(i_node->GetSolutionStepValue(rVariable, SolutionStepNumber));
            float* tmp = (float*)(&data_list[0]);
            vtk_write_compressed ( pResultFile, (char*)tmp, sizeof(*tmp)*data_list.size());
        }

        VTK_fEndResult( pResultFile );

        Timer::Stop("Writing Results");
    }

    /**
     * writes nodal results for variables of type array_1d
     */
    void WriteNodalResults( FILE* pResultFile,
                            Variable<array_1d<double, 3> > const& rVariable,
                            NodesContainerType& rNodes,
                            const std::size_t& SolutionStepNumber ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), 3, mMode );

        if (mMode == VTK_PostAscii)
        {
            for ( NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node)
            {
                array_1d<double, 3>& temp = i_node->GetSolutionStepValue( rVariable, SolutionStepNumber );
                VTK_fWriteVector3( pResultFile, i_node->Id(), temp[0], temp[1], temp[2] );
            }
        }
        else if (mMode == VTK_PostBinary)
        {
            std::vector<float> data_list;
            data_list.reserve(3*(rNodes.end() - rNodes.begin()));
            for ( NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node)
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
     * writes nodal results for variables of type Vector
     */
    void WriteNodalResults( FILE* pResultFile,
                            Variable<Vector> const& rVariable,
                            NodesContainerType& rNodes,
                            const std::size_t& SolutionStepNumber,
                            const std::size_t& NumberOfComponents ) const
    {
        Timer::Start("Writing Results");

        VTK_fBeginResult( pResultFile, (char*)(rVariable.Name().c_str()), NumberOfComponents, mMode );

        if (mMode == VTK_PostAscii)
        {
            double* temp = new double[NumberOfComponents];
            for ( NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node)
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
            for ( NodesContainerType::iterator i_node = rNodes.begin(); i_node != rNodes.end() ; ++i_node)
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

private:
    /**
     * assignment operator
     */
    VtkIO& operator=(VtkIO const& rOther);

    /**
     * Copy constructor
     */
    VtkIO(VtkIO const& rOther);

}; // Class VtkIO

}// namespace Kratos.

#endif // KRATOS_VTK_IO_H_INCLUDED  defined

