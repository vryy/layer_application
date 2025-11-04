/*
    Data Reader and Inspector for GiD binary post results
    Copyright (C) 2016  Hoang-Giang Bui

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//
//   Project Name:        GiD Binary Reader
//   Project Description: Dump the data from GiD binary post-processing file
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: Nov 12, 2012 $
//   Revision:            $Revision: 1.0 $
//   Date:                $Date: Feb 3, 2016 $
//   Revision:            $Revision: 2.0 $
//
//


#if !defined(GIDPOST_READER_H_INCLUDED)
#define GIDPOST_READER_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <vector>
#include <map>

// Project includes
#include "includes/define.h"

namespace Kratos
{

/**
 * Abstract reader for post.xxx result file
*/
class GiDPostReader
{
public:

    /// Type Definitions
    KRATOS_CLASS_POINTER_DEFINITION(GiDPostReader);

    typedef int gid_index_t;
    typedef int gid_size_t;
    typedef float gid_value_t;

    /// Default constructor.
    GiDPostReader(const std::string& ResultsDatafile) : m_filename(ResultsDatafile) {}

    /// Destructor.
    virtual ~GiDPostReader() {}

    /// Get the associated result file name
    std::string GetFileName() const {return m_filename;}

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "GiDPostReader";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}

    /// Get the names of all the meshes in the file
    virtual std::vector<std::string> GetMeshesName() const
    {
        return {};
    }

    /// Get the info of the mesh
    virtual void GetMeshInfo(const std::string& Name, int& Dim, std::string& ElemType) const
    {}

    /// Read the coordinates from mesh name from the binary
    virtual void ReadMesh(const std::string& Name, std::map<int, std::vector<double> >& rCoordinates)
    {}

    /// Read the connectivities from mesh name from the binary
    virtual void ReadMesh(const std::string& Name, std::map<int, std::vector<int> >& rConnectivities)
    {}

    /// Read a specific mesh name from the binary
    void ReadMesh(const std::string& Name, std::map<int, std::vector<double> >& rCoordinates, std::map<int, std::vector<int> >& rConnectivities)
    {
        this->ReadMesh(Name, rCoordinates);
        this->ReadMesh(Name, rConnectivities);
    }

    /// Get the names of all the nodal scalar values
    virtual std::vector<std::string> GetNodalScalarValuesName() const
    {
        return {};
    }

    /// Get the names of all the nodal vector values
    virtual std::vector<std::string> GetNodalVectorValuesName() const
    {
        return {};
    }

    /// Get the names of all the Gauss point scalar values
    virtual std::vector<std::pair<std::string, std::string> > GetGaussPointScalarValuesName() const
    {
        return {};
    }

    /// Get the names of all the Gauss point vector values
    virtual std::vector<std::pair<std::string, std::string> > GetGaussPointVectorValuesName() const
    {
        return {};
    }

    /// Get the names of all the Gauss point matrix values
    virtual std::vector<std::pair<std::string, std::string> > GetGaussPointMatrixValuesName() const
    {
        return {};
    }

    /// Read the nodal values (as scalar)
    /// step_list: all the time steps that simulation produces results
    /// output values: map key is node index
    ///                map values is series of scalar results at node at multiple time step
    virtual void ReadNodalScalarValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<double> >& rValues)
    {}

    /// Read the nodal values (as vector)
    /// step_list: all the time steps that simulation produces results
    /// output values: map key is node index
    ///                map values is series of vector results at node at multiple time step
    virtual void ReadNodalVectorValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues, std::size_t vector_size)
    {}

    // /// Read the nodal values (as vector)
    // /// step_list: all the time steps that simulation produces results
    // /// output values: map key is node index
    // ///                map values is series of vector results at node at multiple time step
    // virtual void ReadNodalVectorValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues, std::size_t vector_size)
    // {}

    virtual void ReadGaussPointRecord(const std::string& GpName)
    {}

    /// Get the Gauss point record info
    virtual void GetGaussPointRecordInfo(const std::string& Name, int& Np, std::string& ElemType) const
    {}

    /// Get the Gauss point record info
    virtual void GetGaussPointRecordInfo(const std::string& GpName, int& Np, std::string& ElemType, std::string& NaturalCoordinates) const
    {}

    /// Get the Gauss point coordinates
    virtual void GetGaussPointRecordCoordinates(const std::string& GpName, std::vector<std::vector<double> >& rCoordinates) const
    {}

    /// Get the Gauss point coordinates
    virtual void GetGaussPointRecordCoordinates(const std::string& GpName, std::vector<array_1d<double, 3> >& rCoordinates) const
    {}

    /// Read the Gauss point values (as scalar)
    /// step_list: all the time steps that simulation produces results
    /// output values: map key is element index
    ///                map values is series of scalar results at Gauss points/element at multiple time step
    virtual void ReadGaussPointScalarValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues)
    {}

    /// Read the Gauss point values (as vector)
    /// step_list: all the time steps that simulation produces results
    /// output values: map key is element index
    ///                map values is series of vector results at Gauss points/element at multiple time step
    virtual void ReadGaussPointVectorValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<std::vector<double> > > >& rValues)
    {}

    /// Read the Gauss point values (as matrix)
    /// step_list: all the time steps that simulation produces results
    /// output values: map key is element index
    ///                map values is series of matrix results at Gauss points/element at multiple time step
    virtual void ReadGaussPointMatrixValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<std::vector<double> > > >& rValues)
    {}

    /// Try to get the directory from file name
    std::string GetDirectory() const
    {
        std::size_t pos = m_filename.find_last_of("/\\"); // handle both '/' and '\'
        std::string directory = (pos == std::string::npos) ? "" : m_filename.substr(0, pos);
        return directory;
    }

    /// Try to get the prefix from file name
    std::string GetPrefix() const
    {
        std::size_t pos = m_filename.rfind('_');
        std::string base = m_filename.substr(0, pos);
        std::size_t pos2 = base.find_last_of("/\\");
        std::string fileName = (pos2 == std::string::npos) ? base : base.substr(pos2 + 1);
        return fileName;
    }

    /// Try to get the time from file name
    double GetTime() const
    {
        std::size_t pos = m_filename.rfind('_');
        std::string last = m_filename.substr(pos + 1);
        std::size_t pos2 = last.find('.');
        std::string last2 = last.substr(pos2 + 1);
        std::size_t pos3 = last2.find('.');
        if (pos3 == std::string::npos)
            return std::stof(last.substr(0, pos2));
        else
            return std::stof(last.substr(0, pos2) + "." + last2.substr(0, pos3));
    }

private:

    std::string m_filename;

}; // Class GiDPostReader

/**
 * output stream function
 */
inline std::ostream& operator << (std::ostream& rOStream, const GiDPostReader& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos.

#endif // GIDPOST_READER_H_INCLUDED  defined
