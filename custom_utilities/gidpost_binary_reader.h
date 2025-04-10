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


#if !defined(GIDPOST_BINARY_READER_H_INCLUDED)
#define GIDPOST_BINARY_READER_H_INCLUDED


// System includes
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <chrono> // c++0x feature

// External includes
#include "boost/algorithm/string.hpp"

#include "zlib.h"
#include "zconf.h"

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "gidpost_reader.h"

namespace Kratos
{

/**
 * Reader for post.bin result file
*/
class GiDPostBinaryReader : public GiDPostReader
{
public:

    /// Type Definitions
    KRATOS_CLASS_POINTER_DEFINITION(GiDPostBinaryReader);

    typedef GiDPostReader BaseType;

    using BaseType::gid_index_t;
    using BaseType::gid_size_t;
    using BaseType::gid_value_t;

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
        bool HasCoord() const {return !(EndCoordPos == StartCoordPos);}
        bool HasElem() const {return !(EndElemPos == StartElemPos);}
        void PrintData(std::ostream& rOStream) const
        {
            rOStream << Name;
            rOStream << " Dim:" << Dim;
            rOStream << " ElemType:" << ElemType;
            rOStream << " Nnode:" << Nnode;
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
        std::string GpName; // name of Gauss point rule, in case the result is at Gauss point
        z_off_t StartPos;
        z_off_t EndPos;
        bool HasResult() const {return !(EndPos == StartPos);}
        void PrintData(std::ostream& rOStream) const
        {
            rOStream << Name;
            rOStream << " Problem:" << Problem;
            rOStream << " Step:" << Step;
            rOStream << " Type:" << Type;
            rOStream << " Location:" << Location;
            rOStream << " GpName:" << GpName;
            rOStream << " StartPos:" << StartPos;
            rOStream << " EndPos:" << EndPos;
        }
    };

    struct GpRecord
    {
        std::string Name;
        std::string ElemType;
        int NumberOfGaussPoints;
        std::string NaturalCoordinates;
        std::vector<std::vector<double> > GpCoordinates;
        z_off_t StartCoordPos;
        z_off_t EndCoordPos;
        void PrintData(std::ostream& rOStream) const
        {
            rOStream << Name;
            rOStream << " ElemType:" << ElemType;
            rOStream << " NumberOfGaussPoints:" << NumberOfGaussPoints;
            rOStream << " NaturalCoordinates:" << NaturalCoordinates;
        }
    };

    GiDPostBinaryReader(const std::string& ResultsDatafile);

    ~GiDPostBinaryReader() override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

    std::vector<std::string> GetMeshesName() const override;

    void GetMeshInfo(const std::string& Name, int& Dim, std::string& ElemType) const override;

    void ReadMesh(const std::string& Name, std::map<int, std::vector<double> >& rCoordinates) override;

    void ReadMesh(const std::string& Name, std::map<int, std::vector<int> >& rConnectivities) override;

    std::vector<std::string> GetNodalScalarValuesName() const override;

    std::vector<std::string> GetNodalVectorValuesName() const override;

    void ReadNodalScalarValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<double> >& rValues) override;

    void ReadNodalVectorValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues, std::size_t vector_size) override;

    /// Get the names of all the Gauss point scalar values
    std::vector<std::pair<std::string, std::string> > GetGaussPointScalarValuesName() const override;

    void ReadGaussPointRecord(const std::string& GpName) override;

    void GetGaussPointRecordInfo(const std::string& GpName, int& Np, std::string& ElemType) const override;

    void GetGaussPointRecordInfo(const std::string& GpName, int& Np, std::string& ElemType, std::string& NaturalCoordinates) const override;

    void GetGaussPointRecordCoordinates(const std::string& GpName, std::vector<std::vector<double> >& rCoordinates) const override;

    void GetGaussPointRecordCoordinates(const std::string& GpName, std::vector<array_1d<double, 3> >& rCoordinates) const override;

    void ReadGaussPointScalarValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues) override;

protected:

    std::vector<MeshHeader> mPosMesh;
    std::vector<ResultHeader> mPosResult;
    std::vector<GpRecord> mGpRecord;

    /* this method extracts strings and parse, could be subjected to errors (be careful!!!) */
    /* MESH shall be the only word and keyword used for mesh */
    /* Result name must come with quotation mark (e.g. Result "DISPLACEMENT"), and no other word Result following this rule but not using for result */
    /* THis is the core sub-routine. It needs to be high performance */
    void ParseAndIndex2();

private:

    gzFile m_file = NULL;

    /* close the file */
    int Close();

    /* open the file to read */
    int Open( const char * name );

    int Open( const std::string name );

    /* initialize reading operation */
    int Init();

    /* Reset the cursor position */
    int Reset();

    /* read a string in GiD binary post file */
    /* this version does not return error code */
//    char* ReadString(char* str);
    void ReadString(char* str);

    int ReadString(char* str, int max_len);

    /* read a string in GiD binary post file */
    /* this version return error code */
    int ReadString(std::string& str);

    int CheckEof();

    z_off_t GetCurrentPosition();

    z_off_t SetCurrentPosition(const z_off_t pos);

    template<class TDataType>
    int Read(TDataType& v)
    {
        return gzread(m_file, &v, sizeof(TDataType));
    }

    void ReadWord(std::string& word);

    /* shift n bytes */
    z_off_t Shift(const z_off_t n);

    bool isStringChar(const char ch);

    bool isCharInString(const char c, const std::string str);

    bool CheckForString(const std::string str);

    bool FindNext(int pos, char&c, const char* str, unsigned int len);

}; // Class GiDPostBinaryReader

}  // namespace Kratos.

#endif // GIDPOST_BINARY_READER_H_INCLUDED  defined
