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


#if !defined(KRATOS_VTK_VTU_IO_H_INCLUDED)
#define  KRATOS_VTK_VTU_IO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"
#include "utilities/timer.h"
#include "custom_io/vtk_io.h"
#include "custom_io/vtk_mesh_container.h"

namespace Kratos
{

/**
 * This class defines an interface to Paraview
 * in order to provide VTK compliant I/O functionality
 * As the name suggests, only the unstructured mesh is supported
 * 14 Aug 2017: supports only nodal scalar & vector results
 */
template<class TMeshContainer = VtkMeshContainer>
class VtkVTUIO : public VtkIO<TMeshContainer>
{
public:
    ///pointer definition of VtkVTUIO
    KRATOS_CLASS_POINTER_DEFINITION(VtkVTUIO);

    ///typedefs
    typedef VtkIO<TMeshContainer> BaseType;
    typedef typename BaseType::MeshType MeshType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;
    typedef typename BaseType::NodesContainerType NodesContainerType;
    typedef typename BaseType::ConditionsContainerType ConditionsContainerType;

    /// Constructor
    VtkVTUIO( const std::string& rDatafilename, const VTK_PostMode& rMode)
    : BaseType( rDatafilename, rMode )
    {
    }

    /// Destructor.
    ~VtkVTUIO() override
    {
    }

    /**
     * This has to be called for each solution step BEFORE any results
     * (on nodes) is written
     * @param name the current solution step (i.e. time)
     * @param rThisMesh the mesh containing the results
     */
    void Initialize( double name, const MeshType& rThisMesh ) override
    {
        if ( !BaseType::mResultFileOpen )
        {
            std::stringstream file_name;
            file_name << BaseType::mResultFileName << std::setprecision(12) << "_" << name << ".vtu";
            BaseType::mResultFile = VTK_fOpenPostVTUFile(file_name.str().c_str(), BaseType::mMode);
            BaseType::mResultFileOpen = true;
        }

        for ( typename std::vector<TMeshContainer>::iterator it = BaseType::MeshContainers().begin();
                        it != BaseType::MeshContainers().end(); it++ )
        {
            it->Reset();
        }

        if ( BaseType::mResultFileOpen )
        {
            for ( typename MeshType::ElementConstantIterator element_iterator = rThisMesh.ElementsBegin();
                    element_iterator != rThisMesh.ElementsEnd(); ++element_iterator)
                for ( typename std::vector<TMeshContainer>::iterator it = BaseType::MeshContainers().begin();
                        it != BaseType::MeshContainers().end(); it++ )
                    if ( it->AddElement( element_iterator ) )
                        break;

            for ( typename MeshType::ConditionConstantIterator conditions_iterator = rThisMesh.ConditionsBegin();
                    conditions_iterator != rThisMesh.ConditionsEnd(); ++conditions_iterator )
                for ( typename std::vector<TMeshContainer>::iterator it = BaseType::MeshContainers().begin();
                        it != BaseType::MeshContainers().end(); it++ )
                    if ( it->AddCondition( conditions_iterator ) )
                        break;

            for ( typename std::vector<TMeshContainer>::iterator it = BaseType::MeshContainers().begin();
                    it != BaseType::MeshContainers().end(); ++it )
            {
                it->FinalizeMeshCreation();
            }
        }
    }

    /**
     * This has to be called for each solution step to write all the results to file
     * have been written
     */
    void Finalize() override
    {
        for ( typename std::vector<TMeshContainer>::iterator it = BaseType::MeshContainers().begin();
                    it != BaseType::MeshContainers().end(); ++it )
        {
            // iterate through element mesh in the mesh container
            for ( typename TMeshContainer::MeshElementsContainerType::const_iterator it2 = it->GetMeshElements().begin();
                    it2 != it->GetMeshElements().end(); ++it2 )
            {
                const ModelPart::ElementsContainerType& MeshElements = it2->second;
                const ModelPart::NodesContainerType& MeshNodes = it->GetMeshElementNodes(it2->first);
                const std::string& MeshName = it->GetMeshElementsName(it2->first);

                // start writing the piece
                VTK_fBeginMesh( BaseType::mResultFile, MeshName.c_str(), MeshNodes.size(), MeshElements.size() );

                // write the points & connectivities
                it->WriteMesh( BaseType::mResultFile, MeshNodes, MeshElements, false, BaseType::mMode );

                // write the nodal results
                BaseType::BeginNodalResultsHeader( BaseType::mResultFile );

                BaseType::WriteNodalResults( BaseType::mResultFile, MeshNodes );

                BaseType::EndNodalResultsHeader( BaseType::mResultFile );

                // write the cell results
                BaseType::BeginCellResultsHeader( BaseType::mResultFile );

                BaseType::WriteCellResults( BaseType::mResultFile, MeshElements );

                BaseType::EndCellResultsHeader( BaseType::mResultFile );

                // end writing the piece
                VTK_fEndMesh ( BaseType::mResultFile );
            }

            // iterate through element mesh in the mesh container
            for ( typename TMeshContainer::MeshConditionsContainerType::const_iterator it2 = it->GetMeshConditions().begin();
                    it2 != it->GetMeshConditions().end(); ++it2 )
            {
                const ModelPart::ConditionsContainerType& MeshConditions = it2->second;
                const ModelPart::NodesContainerType& MeshNodes = it->GetMeshConditionNodes(it2->first);
                const std::string& MeshName = it->GetMeshConditionsName(it2->first);

                // start writing the piece
                VTK_fBeginMesh( BaseType::mResultFile, MeshName.c_str(), MeshNodes.size(), MeshConditions.size() );

                // write the points & connectivities
                it->WriteMesh( BaseType::mResultFile, MeshNodes, MeshConditions, false, BaseType::mMode );

                // write the nodal results
                BaseType::BeginNodalResultsHeader( BaseType::mResultFile );

                BaseType::WriteNodalResults( BaseType::mResultFile, MeshNodes );

                BaseType::EndNodalResultsHeader( BaseType::mResultFile );

                // write the cell results
                BaseType::BeginCellResultsHeader( BaseType::mResultFile );

                BaseType::WriteCellResults( BaseType::mResultFile, MeshConditions );

                BaseType::EndCellResultsHeader( BaseType::mResultFile );

                // end writing the piece
                VTK_fEndMesh ( BaseType::mResultFile );
            }
        }

        // close the result file
        VTK_fClosePostVTUFile( BaseType::mResultFile );
        BaseType::mResultFileOpen = false;
    }

private:
    /**
     * assignment operator
     */
    VtkVTUIO& operator=(VtkVTUIO const& rOther);

    /**
     * Copy constructor
     */
    VtkVTUIO(VtkVTUIO const& rOther);

}; // Class VtkVTUIO

} // namespace Kratos.

#endif // KRATOS_VTK_VTU_IO_H_INCLUDED  defined
