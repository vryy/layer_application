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


#if !defined(KRATOS_VTK_VTM_IO_H_INCLUDED)
#define  KRATOS_VTK_VTM_IO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_io/vtk_io.h"
#include "custom_io/vtk_mesh_container.h"

namespace Kratos
{

/**
 * This class defines an interface to Paraview
 * in order to provide VTK compliant I/O functionality
 * 14 Aug 2017: supports only nodal scalar & vector results
 * 26 Jan 2020: this class supports vtkMultiBlockDataSet in vtm format
 */
template<class TMeshContainer = VtkMeshContainer>
class VtkVTMIO : public VtkIO<TMeshContainer>
{
public:
    ///pointer definition of VtkVTMIO
    KRATOS_CLASS_POINTER_DEFINITION(VtkVTMIO);

    ///typedefs
    typedef VtkIO<TMeshContainer> BaseType;
    typedef typename BaseType::MeshType MeshType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;
    typedef typename BaseType::NodesContainerType NodesContainerType;
    typedef typename BaseType::ConditionsContainerType ConditionsContainerType;

    /// Constructor
    VtkVTMIO( const std::string& rDatafilename, const VTK_PostMode& rMode )
    : BaseType( rDatafilename, rMode )
    {
    }

    /// Destructor.
    ~VtkVTMIO() override
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
        mName = name;

        if ( !BaseType::mResultFileOpen )
        {
            std::stringstream file_name;
            file_name << BaseType::mResultFileName << std::setprecision(12) << "_" << name << ".vtm";
            BaseType::mResultFile = VTK_fOpenPostVTMFile(file_name.str().c_str(), BaseType::mMode);
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
        FILE* tmp_file;

        std::string base_filename = BaseType::mResultFileName.substr(BaseType::mResultFileName.find_last_of("/\\") + 1);

        unsigned int file_index = 0;
        for ( typename std::vector<TMeshContainer>::iterator it = BaseType::MeshContainers().begin();
                    it != BaseType::MeshContainers().end(); ++it )
        {
            // iterate through element mesh in the mesh container
            for ( typename TMeshContainer::MeshElementsContainerType::const_iterator it2 = it->GetMeshElements().begin();
                    it2 != it->GetMeshElements().end(); ++it2 )
            {
                const std::string& MeshName = it->GetMeshElementsName(it2->first);

                std::stringstream tmp_filename;
                tmp_filename << BaseType::mResultFileName << "_" << MeshName
                              << std::setprecision(12) << "_" << mName << ".vtu";
                tmp_file = VTK_fOpenPostVTUFile(tmp_filename.str().c_str(), BaseType::mMode);

                const ModelPart::ElementsContainerType& MeshElements = it2->second;
                const ModelPart::NodesContainerType& MeshNodes = it->GetMeshElementNodes(it2->first);

                // start writing the piece
                VTK_fBeginMesh( tmp_file, MeshName.c_str(), MeshNodes.size(), MeshElements.size() );

                // write the points & connectivities
                it->WriteMesh( tmp_file, MeshNodes, MeshElements, false, BaseType::mMode );

                // write the nodal results
                BaseType::BeginNodalResultsHeader( tmp_file );

                BaseType::WriteNodalResults( tmp_file, MeshNodes );

                BaseType::EndNodalResultsHeader( tmp_file );

                // write the cell results
                BaseType::BeginCellResultsHeader( tmp_file );

                BaseType::WriteCellResults( tmp_file, MeshElements );

                BaseType::EndCellResultsHeader( tmp_file );

                // end writing the piece
                VTK_fEndMesh ( tmp_file );

                // close the vtu file
                VTK_fClosePostVTUFile( tmp_file );

                // register the entry to vtm file
                std::stringstream tmp_filename_short;
                tmp_filename_short << base_filename << "_" << MeshName
                                   << std::setprecision(12) << "_" << mName << ".vtu";
                VTK_fBeginDataset( BaseType::mResultFile, MeshName.c_str(), tmp_filename_short.str().c_str(), file_index++ );
                VTK_fEndDataset( BaseType::mResultFile );
            }

            // iterate through condition mesh in the mesh container
            for ( typename TMeshContainer::MeshConditionsContainerType::const_iterator it2 = it->GetMeshConditions().begin();
                    it2 != it->GetMeshConditions().end(); ++it2 )
            {
                const std::string& MeshName = it->GetMeshConditionsName(it2->first);

                std::stringstream tmp_filename;
                tmp_filename << BaseType::mResultFileName << "_" << MeshName
                              << std::setprecision(12) << "_" << mName << ".vtu";
                tmp_file = VTK_fOpenPostVTUFile(tmp_filename.str().c_str(), BaseType::mMode);

                const ModelPart::ConditionsContainerType& MeshConditions = it2->second;
                const ModelPart::NodesContainerType& MeshNodes = it->GetMeshConditionNodes(it2->first);

                // start writing the piece
                VTK_fBeginMesh( tmp_file, MeshName.c_str(), MeshNodes.size(), MeshConditions.size() );

                // write the points & connectivities
                it->WriteMesh( tmp_file, MeshNodes, MeshConditions, false, BaseType::mMode );

                // write the nodal results
                BaseType::BeginNodalResultsHeader( tmp_file );

                BaseType::WriteNodalResults( tmp_file, MeshNodes );

                BaseType::EndNodalResultsHeader( tmp_file );

                // write the cell results
                BaseType::BeginCellResultsHeader( tmp_file );

                BaseType::WriteCellResults( tmp_file, MeshConditions );

                BaseType::EndCellResultsHeader( tmp_file );

                // end writing the piece
                VTK_fEndMesh ( tmp_file );

                // close the vtu file
                VTK_fClosePostVTUFile( tmp_file );

                // register the entry to vtm file
                std::stringstream tmp_filename_short;
                tmp_filename_short << base_filename << "_" << MeshName
                                   << std::setprecision(12) << "_" << mName << ".vtu";
                VTK_fBeginDataset( BaseType::mResultFile, MeshName.c_str(), tmp_filename_short.str().c_str(), file_index++ );
                VTK_fEndDataset( BaseType::mResultFile );
            }
        }

        // close the result file
        VTK_fClosePostVTMFile( BaseType::mResultFile );
        BaseType::mResultFileOpen = false;
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const override
    {
        return "VTK VTM IO";
    }

private:

    double mName;

    /**
     * assignment operator
     */
    VtkVTMIO& operator=(VtkVTMIO const& rOther);

    /**
     * Copy constructor
     */
    VtkVTMIO(VtkVTMIO const& rOther);

}; // Class VtkVTMIO

} // namespace Kratos.

#endif // KRATOS_VTK_VTM_IO_H_INCLUDED  defined
