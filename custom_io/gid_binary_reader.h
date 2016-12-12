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


#if !defined(GID_BINARY_READER_H_INCLUDED)
#define GID_BINARY_READER_H_INCLUDED


// System includes
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <chrono> // c++0x feature



// External includes
#include "boost/smart_ptr.hpp"
#include "boost/algorithm/string.hpp"

#include "zlib.h"
#include "zconf.h"

#define GZ_WATCH(a) if(a <= 0) return a;
#define WATCH(a) std::cout << #a << ": " << a << std::endl;
#define ERROR(a, b) {std::cerr << a << " " << b << std::endl; exit(0);}

namespace Kratos
{

/**
Detail class definition.
*/
class GiDBinaryReader
{
public:

    /// Type Definitions
    typedef boost::shared_ptr<GiDBinaryReader> Pointer;

    typedef int gid_index_t;
    typedef int gid_size_t;
    typedef float gid_value_t;

    struct MeshHeader
    {
        std::string Name;
        int Dim;
        std::string ElemType;
        int Nnode;
        z_off_t StartCoordPos;
        z_off_t EndCoordPos;
        z_off_t StartElemPos;
        z_off_t EndElemPos;
        void PrintData(std::ostream& rOStream) const
        {
            rOStream << Name;
            rOStream << " Dim:" << Dim;
            rOStream << " " << ElemType;
            rOStream << " " << Nnode;
            rOStream << " StartCoordPos:" << StartCoordPos;
            rOStream << " EndCoordPos:" << EndCoordPos;
            rOStream << " StartElemPos:" << StartElemPos;
            rOStream << " EndElemPos:" << EndElemPos;
        }
    };

    struct ResultHeader
    {
        std::string Name;
        std::string Problem;
        double Step;
        std::string Type;
        std::string Location;
        z_off_t StartPos;
        z_off_t EndPos;
        void PrintData(std::ostream& rOStream) const
        {
            rOStream << Name;
            rOStream << " " << Problem;
            rOStream << " Step:" << Step;
            rOStream << " " << Type;
            rOStream << " " << Location;
            rOStream << " StartPos:" << StartPos;
            rOStream << " EndPos:" << EndPos;
        }
    };

    /// Default constructor.
    GiDBinaryReader(const std::string& ResultsDatafile) : m_filename(ResultsDatafile)
    {
        int fail = Open(m_filename.c_str());

        if(fail != 0)
            ERROR("Fail to open file", m_filename);

        Init();
        ParseAndIndex2();
//        PrintData(std::cout);
    }

    /// Destructor.
    virtual ~GiDBinaryReader()
    {
        int fail = Close();
        if(fail == Z_OK)
        {
            std::cout << "File is closed successfully" << std::endl;
            std::cout << "-------------------------------------------------------------------" << std::endl;
            std::cout << "|                           GLÜCK WÜNSCH                          |" << std::endl;
            std::cout << "-------------------------------------------------------------------" << std::endl;
        }
        else
            std::cout << "Error closing file: " << fail << std::endl;
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Data Reader and Inspector for GiD binary post results, (c) Hoang-Giang Bui 2016";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "There is(are) " << mPosMesh.size() << " mesh(es) in the binary data" << std::endl;
        rOStream << "There is(are) " << mPosResult.size() << " result(s) in the binary data" << std::endl;
        rOStream << "Mesh summary:" << std::endl;
        for(std::size_t i = 0; i < mPosMesh.size(); ++i)
        {
            mPosMesh[i].PrintData(rOStream << "  ");
            rOStream << std::endl;
        }
        rOStream << "Result summary:" << std::endl;
        for(std::size_t i = 0; i < mPosResult.size(); ++i)
        {
            mPosResult[i].PrintData(rOStream << "  ");
            rOStream << std::endl;
        }
    }

    /// Get the names of all the meshes in the file
    std::vector<std::string> GetMeshesName() const
    {
        std::vector<std::string> mesh_list;
        for(std::size_t i = 0; i < mPosMesh.size(); ++i)
        {
            mesh_list.push_back(mPosMesh[i].Name);
        }
        return mesh_list;
    }

    /// step_list: all the time steps that simulation produces results
    /// output values: map key is node index
    ///                map values is series of scalar results at node at multiple time step
    void ReadNodalScalarValues(std::string Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<double> >& rValues)
    {
        Reset();

        step_list.clear();
        rValues.clear();

        for(std::size_t i = 0; i < mPosResult.size(); ++i)
        {
            std::string RefName = '"' + Name + '"';

            bool check = (mPosResult[i].Name.compare(RefName) == 0) || (mPosResult[i].Name.compare(Name) == 0);
            check = check && (mPosResult[i].Type.compare("Scalar") == 0);
            check = check && (mPosResult[i].Location.compare("OnNodes") == 0);

            if(check)
            {
                step_list.push_back(mPosResult[i].Step);

                SetCurrentPosition(mPosResult[i].StartPos);

                std::string line;
                ReadString(line);

                std::vector<std::string> fields;
                boost::split(fields, line, boost::is_any_of(" \n\0"));
//                std::cout << "fields:";
//                for(std::size_t j = 0; j < fields.size(); ++j)
//                    std::cout << " " << fields[j];
//                std::cout << std::endl;

                if(fields[0].compare("Values") == 0)
                {
                    z_off_t indexed = atoi(fields[1].c_str());
                    gid_index_t id;
                    gid_value_t v;
                    z_off_t CurPos;
//                    WATCH(GetCurrentPosition())

                    if(indexed == -1)
                    {
                        do
                        {
                            Read(id);
                            Read(v);
                            if(id > 0)
                                rValues[id].push_back(v);
                            // REMARKS: for some reason, the last line contains the value pair: [id, v] = [-1, ##]. Probably to terminate the sequence of value. I don't know how to exclude this when I compute the end cursor position, so I make a naive comparison to make sure correct values are filtered out.
                            CurPos = GetCurrentPosition();
                        }
                        while(CurPos < mPosResult[i].EndPos);
                    }
                    else
                    {
                        for(z_off_t i = 0; i < indexed; ++i)
                        {
                            Read(id);
                            Read(v);
                            rValues[id].push_back(v);
                        }
                    }
                }
            }
        }
    }

    /// TODO
    void ReadNodalVectorValues(std::string Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues, std::size_t vector_size)
    {
        Reset();

        step_list.clear();
        rValues.clear();

        for(std::size_t i = 0; i < mPosResult.size(); ++i)
        {
            std::string RefName = '"' + Name + '"';

            bool check = (mPosResult[i].Name.compare(RefName) == 0 || mPosResult[i].Name.compare(Name) == 0);
            check = check && (mPosResult[i].Type.compare("Vector") == 0);
            check = check && (mPosResult[i].Location.compare("OnNodes") == 0);

            if(check)
            {
                step_list.push_back(mPosResult[i].Step);

                SetCurrentPosition(mPosResult[i].StartPos);

                std::string line;
                ReadString(line);

                std::vector<std::string> fields;
                boost::split(fields, line, boost::is_any_of(" \n\0"));
//                std::cout << "fields:";
//                for(std::size_t j = 0; j < fields.size(); ++j)
//                    std::cout << " " << fields[j];
//                std::cout << std::endl;

                if(fields[0].compare("Values") == 0)
                {
                    z_off_t last_node_id = atoi(fields[1].c_str());
                    gid_index_t id;
                    z_off_t CurPos;
//                    WATCH(GetCurrentPosition())

                    if(last_node_id == -1)
                    {
                        do
                        {
                            Read(id);
                            std::vector<double> v(vector_size);
                            for(std::size_t j = 0; j < vector_size; ++j)
                            {
                                gid_value_t t;
                                Read(t);
                                v[j] = t;
                            }
                            if(id > 0)
                                rValues[id].push_back(v);
                            // REMARKS: for some reason, the last line contains the value pair: [id, v] = [-1, ##]. Probably to terminate the sequence of value. I don't know how to exclude this when I approach the end cursor position, so I make a naive comparison to make sure correct values are filtered out.
                            CurPos = GetCurrentPosition();
                        }
                        while(CurPos < mPosResult[i].EndPos);
                    }
                    else
                    {
                        for(z_off_t i = 0; i < last_node_id; ++i)
                        {
                            Read(id);
                            std::vector<double> v(vector_size);
                            for(std::size_t j = 0; j < vector_size; ++j)
                            {
                                gid_value_t t;
                                Read(t);
                                v[j] = t;
                            }
                            rValues[id].push_back(v);
                        }
                    }
                }
            }
        }
    }

    /// Read a specific mesh name from the binary
    void ReadMesh(std::string Name, std::map<int, std::vector<double> >& rCoordinates,
                  std::map<int, std::vector<int> >& rConnectivities)
    {
        Reset();

        std::string RefName = '"' + Name + '"';
        for(std::size_t i = 0; i < mPosMesh.size(); ++i)
        {
            bool check = (mPosMesh[i].Name.compare(RefName) == 0
                    || mPosMesh[i].Name.compare(Name) == 0);

            if(check)
            {
                bool coordinates_are_read = false;
                bool connectivities_are_read = false;
                std::string line;
                std::vector<std::string> fields;
                z_off_t CurPos;

                /* iterate through all the coordinates */
                SetCurrentPosition(mPosMesh[i].StartCoordPos);
                do
                {
                    ReadString(line);
//                    WATCH(line)

                    fields.clear();
                    boost::split(fields, line, boost::is_any_of(" "));

                    if(fields[0].compare("Coordinates") == 0)
                    {
                        z_off_t last_node_id = atoi(fields[1].c_str());
                        gid_index_t id;
                        gid_value_t x, y, z;

                        if(last_node_id == -1)
                        {
                            while(true)
                            {
                                Read(id);
                                CurPos = GetCurrentPosition();
                                if(CurPos >= mPosMesh[i].EndCoordPos)
                                {
                                    break;
                                }
                                std::vector<double> v(3);
                                for(std::size_t j = 0; j < 3; ++j)
                                {
                                    gid_value_t t;
                                    Read(t);
                                    v[j] = t;
                                }
                                if(id > 0)
                                    rCoordinates[id] = v;
                            }
                        }
                        else
                        {
                            for(z_off_t k = 0; k < last_node_id; ++k)
                            {
                                Read(id);
                                std::vector<double> v(3);
                                for(std::size_t j = 0; j < 3; ++j)
                                {
                                    gid_value_t t;
                                    Read(t);
                                    v[j] = t;
                                }
                                rCoordinates[id] = v;
                            }
                        }
                        coordinates_are_read = true;
                    }
                }
                while(!(coordinates_are_read));

                /* iterate through all the connectivities */
                SetCurrentPosition(mPosMesh[i].StartElemPos);
                do
                {
                    ReadString(line);
//                    WATCH(line)

                    fields.clear();
                    boost::split(fields, line, boost::is_any_of(" "));

                    if(fields[0].compare("Elements") == 0)
                    {
                        z_off_t last_elem_id = atoi(fields[1].c_str());
                        gid_index_t id;
                        gid_value_t x, y, z;

                        if(last_elem_id == -1)
                        {
                            while(true)
                            {
                                Read(id);
                                CurPos = GetCurrentPosition();
                                if(CurPos >= mPosMesh[i].EndElemPos)
                                {
                                    break;
                                }
                                std::vector<gid_index_t> e(mPosMesh[i].Nnode);
                                for(std::size_t j = 0; j < e.size(); ++j)
                                {
                                    Read(e[j]);
                                }
                                gid_index_t dummy; // TODO: check what is this
                                Read(dummy);
                                if(id > 0)
                                    rConnectivities[id] = e;
                            }
                        }
                        else
                        {
                            for(z_off_t k = 0; k < last_elem_id; ++k)
                            {
                                Read(id);
                                std::vector<gid_index_t> e(mPosMesh[i].Nnode);
                                for(std::size_t j = 0; j < e.size(); ++j)
                                {
                                    Read(e[j]);
                                }
                                gid_index_t dummy; // TODO: check what is this
                                Read(dummy);
                                if(id > 0)
                                    rConnectivities[id] = e;
                            }
                        }
                        connectivities_are_read = true;
                    }
                }
                while(!(connectivities_are_read));
            }
        }
    }

protected:

    std::vector<MeshHeader> mPosMesh;
    std::vector<ResultHeader> mPosResult;

    /* this method extracts strings and parse, could be subjected to errors (be careful!!!) */
    /* MESH shall be the only word and keyword used for mesh */
    /* Result name must come with quotation mark (e.g. Result "DISPLACEMENT"), and no other word Result following this rule but not using for result */
    /* THis is the core sub-routine. It needs to be high performance */
    void ParseAndIndex2()
    {
        auto time_begin = std::chrono::high_resolution_clock::now();

        char c;
        std::map<std::size_t, std::vector<z_off_t> > index;

        std::vector<std::string> keywords;
        keywords.push_back("MESH");
        keywords.push_back("Result");

        // search for the keywords
        while(!CheckEof())
        {
            Read(c);

            for(std::size_t i = 0; i < keywords.size(); ++i)
            {
                std::size_t cnt = 1;
                bool detect = false;
                while(c == keywords[i][cnt-1])
                {
                    if(cnt == keywords[i].size())
                    {
                        detect = true;
                        index[i].push_back(GetCurrentPosition() - cnt);
                        break;
                    }
                    Read(c);
                    ++cnt;
                }
            }
        }

        // handle keyword "MESH"
        for(std::size_t i = 0; i < index[0].size(); ++i)
        {
            // set the cursor position
            SetCurrentPosition(index[0][i] + 5);

            // read the header
            MeshHeader mesh;
            ReadWord(mesh.Name);

            std::string dimension;
            ReadWord(dimension);

            std::string tmp;
            ReadWord(tmp);
            mesh.Dim = atoi(tmp.c_str());

            std::string elem_type;
            ReadWord(elem_type);
            ReadWord(mesh.ElemType);

            std::string nnode;
            ReadWord(nnode);
            tmp.clear();
            ReadWord(tmp);
            mesh.Nnode = atoi(tmp.c_str());

            mesh.StartCoordPos = GetCurrentPosition();

            // search for the next End keyword
            bool found = false;
            while(!CheckEof())
            {
                Read(c);
                if(c == 'E' || c == 'e')
                {
                    Read(c);
                    if(c == 'n' || c == 'N')
                    {
                        Read(c);
                        if(c == 'd' || c == 'D')
                        {
                            Read(c);
                            if(c == ' ')
                            {
                                Read(c);
                                if(c == 'c' || c == 'C')
                                {
                                    Read(c);
                                    if(c == 'o' || c == 'O')
                                    {
                                        Read(c);
                                        if(c == 'o' || c == 'O')
                                        {
                                            Read(c);
                                            if(c == 'r' || c == 'R')
                                            {
                                                Read(c);
                                                if(c == 'd' || c == 'D')
                                                {
                                                    Read(c);
                                                    if(c == 'i' || c == 'I')
                                                    {
                                                        Read(c);
                                                        if(c == 'n' || c == 'N')
                                                        {
                                                            Read(c);
                                                            if(c == 'a' || c == 'A')
                                                            {
                                                                Read(c);
                                                                if(c == 't' || c == 'T')
                                                                {
                                                                    Read(c);
                                                                    if(c == 'e' || c == 'E')
                                                                    {
                                                                        Read(c);
                                                                        if(c == 's' || c == 'S')
                                                                        {
                                                                            found = true;
                                                                            break;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                         }
                                     }
                                 }
                            }
                        }
                    }
                }
            }
            if(found)
                mesh.EndCoordPos = GetCurrentPosition() - 15*sizeof(char);
            else
                continue;
//                ERROR("That should not happen", "No matching End Coordinates keyword")
            // read through the coordinates
            std::string dummy;
            ReadWord(dummy);

            mesh.StartElemPos = GetCurrentPosition();
            // search for the next End keyword
            found = false;
            while(!CheckEof())
            {
                Read(c);
                if(c == 'E' || c == 'e')
                {
                    Read(c);
                    if(c == 'n' || c == 'N')
                    {
                        Read(c);
                        if(c == 'd' || c == 'D')
                        {
                            Read(c);
                            if(c == ' ')
                            {
                                Read(c);
                                if(c == 'e' || c == 'E')
                                {
                                    Read(c);
                                    if(c == 'l' || c == 'L')
                                    {
                                        Read(c);
                                        if(c == 'e' || c == 'E')
                                        {
                                            Read(c);
                                            if(c == 'm' || c == 'M')
                                            {
                                                Read(c);
                                                if(c == 'e' || c == 'E')
                                                {
                                                    Read(c);
                                                    if(c == 'n' || c == 'N')
                                                    {
                                                        Read(c);
                                                        if(c == 't' || c == 'T')
                                                        {
                                                            Read(c);
                                                            if(c == 's' || c == 'S')
                                                            {
                                                                found = true;
                                                                break;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(found)
                mesh.EndElemPos = GetCurrentPosition() - 12*sizeof(char);
            else
                continue;
//                ERROR("That should not happen", "No matching End Elements keyword")

            mPosMesh.push_back(mesh);
        }

        // handle keyword "Result"
        for(std::size_t i = 0; i < index[1].size(); ++i)
        {
            // set the cursor position
            SetCurrentPosition(index[1][i] + 7);

            // read the header
            ResultHeader result;

            ReadWord(result.Name);
            ReadWord(result.Problem);

            std::string step_str;
            ReadWord(step_str);
            result.Step = atof(step_str.c_str());

            ReadWord(result.Type);
            ReadWord(result.Location);

            result.StartPos = GetCurrentPosition();

            // search for the next End keyword
            bool found = false;
            while(!CheckEof())
            {
                Read(c);
                if(c == 'E' || c == 'e')
                {
                    Read(c);
                    if(c == 'n' || c == 'N')
                    {
                        Read(c);
                        if(c == 'd' || c == 'D')
                        {
                            Read(c);
                            if(c == ' ')
                            {
                                Read(c);
                                if(c == 'v' || c == 'V')
                                {
                                    Read(c);
                                    if(c == 'a' || c == 'A')
                                    {
                                        Read(c);
                                        if(c == 'l' || c == 'L')
                                        {
                                            Read(c);
                                            if(c == 'u' || c == 'U')
                                            {
                                                Read(c);
                                                if(c == 'e' || c == 'E')
                                                {
                                                    Read(c);
                                                    if(c == 's' || c == 'S')
                                                    {
                                                        found = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(found)
                result.EndPos = GetCurrentPosition() - 10*sizeof(char);
            else
                continue;
//                ERROR("That should not happen", "No matching End keyword")

            mPosResult.push_back(result);
        }

        auto time_end = std::chrono::high_resolution_clock::now();
        auto elapsed = time_end - time_begin;
        long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        std::cout << "Indexing " + m_filename + " completed... " << microseconds << "us" << std::endl;
    }

private:

    gzFile m_file = NULL;

    /* close the file */
    int Close()
    {
        int fail = Z_OK;
        if(m_file)
        {
            fail = gzclose(m_file);
            m_file = NULL;
        }
        return fail;
    }

    /* open the file to read */
    int Open( const char * name )
    {
        Close();
        m_file = gzopen(name, "rb");
        return m_file == NULL;
    }

    int Open( const std::string name )
    {
        return Open(name.c_str());
    }

    /* initialize reading operation */
    int Init()
    {
        /* skip the first int */
        return gzseek(m_file, sizeof(int), SEEK_SET);
    }

    /* Reset the cursor position */
    int Reset()
    {
        return gzrewind(m_file);
    }

    /* read a string in GiD binary post file */
    /* this version does not return error code */
    char* ReadString(char* str)
    {
        int len;
        gzread(m_file, &len, sizeof(int));
        if(strlen(str) < len)
            str = (char *)malloc(sizeof(char) * len);
        gzread(m_file, str, len*sizeof(char));
        return str;
    }

    /* read a string in GiD binary post file */
    /* this version return error code */
    int ReadString(std::string& str)
    {
        int len;
        int error_code = gzread(m_file, &len, sizeof(int)); // read the length of the string
        GZ_WATCH(error_code);
        if(len > 0)
        {
            if(str.size() != len)
                str.resize(len);
            error_code = gzread(m_file, (char *)str.c_str(), len*sizeof(char));
            GZ_WATCH(error_code);
            return 0;
        }
        else
        {
            str.resize(0);
            return 1;
        }
    }

    int CheckEof()
    {
        return gzeof(m_file);
    }

    z_off_t GetCurrentPosition()
    {
        return gztell(m_file);
    }

    z_off_t SetCurrentPosition(const z_off_t pos)
    {
        return gzseek(m_file, pos, SEEK_SET);
    }

    template<class TDataType>
    int Read(TDataType& v)
    {
        return gzread(m_file, &v, sizeof(TDataType));
    }

    int ReadWord(std::string& word)
    {
        char c;
        Read(c);
        while((c != ' ') && (c != '\0'))
        {
            word.push_back(c);
            Read(c);
        }
    }

    /* shift n bytes */
    z_off_t Shift(const z_off_t n)
    {
        return gzseek(m_file, n, SEEK_CUR);
    }

    bool isStringChar(const char ch) {
        if (ch >= 'a' && ch <= 'z')
            return true;
        if (ch >= 'A' && ch <= 'Z')
            return true;
        return false;
    }

    bool isCharInString(const char c, const std::string str)
    {
        for(int i = 0; i < str.size(); ++i)
        {
            if(c == str[i])
                return true;
        }
        return false;
    }

    bool CheckForString(const std::string str)
    {
        char c;
        for(int i = 0; i < str.size(); ++i)
        {
            Read(c);
            if(c != str[i])
            {
                Shift(-(i+1)*sizeof(char));
                return false;
            }
        }
        Shift(-str.size()*sizeof(char));
        return true;
    }

    std::string m_filename;

}; // Class GiDBinaryReader

/**
 * output stream function
 */
inline std::ostream& operator << (std::ostream& rOStream, const GiDBinaryReader& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}


}  // namespace Kratos.

#endif // GID_BINARY_READER_H_INCLUDED  defined 

