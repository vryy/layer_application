#include <set>

#include "gidpost_binary_reader.h"

#define GZ_WATCH(a) if(a <= 0) return a;
#define WATCH(a) std::cout << #a << ": " << a << std::endl;
#define ERROR(a, b) {std::cerr << a << " " << b << std::endl; exit(0);}
#define ERROR2(a, b, c) {std::cerr << a << " " << b << " " << c << std::endl; exit(0);}

namespace Kratos
{

GiDPostBinaryReader::GiDPostBinaryReader(const std::string& ResultsDatafile) : BaseType(ResultsDatafile)
{
    std::cout << "GiDPostBinaryReader, Zlib version = " << ZLIB_VERSION << std::endl;

    int fail = Open(BaseType::GetFileName().c_str());

    if(fail != 0)
        ERROR("Fail to open file", BaseType::GetFileName());

    Init();
    ParseAndIndex2();
//        PrintData(std::cout);
}

GiDPostBinaryReader::~GiDPostBinaryReader()
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

std::string GiDPostBinaryReader::Info() const
{
    return "Data Reader and Inspector for GiD binary post results, (c) Hoang-Giang Bui 2016-2022";
}

void GiDPostBinaryReader::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void GiDPostBinaryReader::PrintData(std::ostream& rOStream) const
{
    rOStream << "There is(are) " << mPosMesh.size() << " mesh(es) in the binary data" << std::endl;
    rOStream << "There is(are) " << mPosResult.size() << " result(s) in the binary data" << std::endl;
    rOStream << "MESH summary:" << std::endl;
    for(std::size_t i = 0; i < mPosMesh.size(); ++i)
    {
        mPosMesh[i].PrintData(rOStream << "  ");
        rOStream << std::endl;
    }
    rOStream << "GaussPoints summary:" << std::endl;
    for(std::size_t i = 0; i < mGpRecord.size(); ++i)
    {
        mGpRecord[i].PrintData(rOStream << "  ");
        rOStream << std::endl;
    }
    rOStream << "Result summary:" << std::endl;
    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        mPosResult[i].PrintData(rOStream << "  ");
        rOStream << std::endl;
    }
}

std::vector<std::string> GiDPostBinaryReader::GetMeshesName() const
{
    std::vector<std::string> mesh_list;
    for(std::size_t i = 0; i < mPosMesh.size(); ++i)
    {
        mesh_list.push_back(mPosMesh[i].Name);
    }
    return mesh_list;
}

void GiDPostBinaryReader::GetMeshInfo(const std::string& Name, int& Dim, std::string& ElemType) const
{
    for(std::size_t i = 0; i < mPosMesh.size(); ++i)
    {
        if (mPosMesh[i].Name.compare(Name) == 0)
        {
            Dim = mPosMesh[i].Dim;
            ElemType = mPosMesh[i].ElemType;
            return;
        }
    }

    KRATOS_ERROR << "Mesh " << Name << " is not found";
}

std::vector<std::string> GiDPostBinaryReader::GetNodalScalarValuesName() const
{
    std::set<std::string> result_list;
    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        if (mPosResult[i].Type.compare("Scalar") == 0
            && (mPosResult[i].Location.compare("OnNodes") == 0))
        {
            result_list.insert(mPosResult[i].Name);
        }
    }
    return std::vector<std::string>(result_list.begin(), result_list.end());
}

std::vector<std::string> GiDPostBinaryReader::GetNodalVectorValuesName() const
{
    std::set<std::string> result_list;
    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        if (mPosResult[i].Type.compare("Vector") == 0
            && (mPosResult[i].Location.compare("OnNodes") == 0))
        {
            result_list.insert(mPosResult[i].Name);
        }
    }
    return std::vector<std::string>(result_list.begin(), result_list.end());
}

void GiDPostBinaryReader::ReadNodalScalarValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<double> >& rValues)
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

void GiDPostBinaryReader::ReadNodalVectorValues(const std::string& Name, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues, std::size_t vector_size)
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
                gid_value_t t;
//                    WATCH(GetCurrentPosition())

                if(last_node_id == -1)
                {
                    do
                    {
                        Read(id);
                        std::vector<double> v(vector_size);
                        for(std::size_t j = 0; j < vector_size; ++j)
                        {
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

std::vector<std::pair<std::string, std::string> > GiDPostBinaryReader::GetGaussPointScalarValuesName() const
{
    std::vector<std::pair<std::string, std::string> > result_list;
    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        if (mPosResult[i].Type.compare("Scalar") == 0
            && (mPosResult[i].Location.compare("OnGaussPoints") == 0))
        {
            result_list.push_back(std::make_pair(mPosResult[i].Name, mPosResult[i].GpName));
        }
    }
    return result_list;
}

std::vector<std::pair<std::string, std::string> > GiDPostBinaryReader::GetGaussPointVectorValuesName() const
{
    std::vector<std::pair<std::string, std::string> > result_list;
    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        if (mPosResult[i].Type.compare("Vector") == 0
            && (mPosResult[i].Location.compare("OnGaussPoints") == 0))
        {
            result_list.push_back(std::make_pair(mPosResult[i].Name, mPosResult[i].GpName));
        }
    }
    return result_list;
}

std::vector<std::pair<std::string, std::string> > GiDPostBinaryReader::GetGaussPointMatrixValuesName() const
{
    std::vector<std::pair<std::string, std::string> > result_list;
    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        if (mPosResult[i].Type.compare("Matrix") == 0
            && (mPosResult[i].Location.compare("OnGaussPoints") == 0))
        {
            result_list.push_back(std::make_pair(mPosResult[i].Name, mPosResult[i].GpName));
        }
    }
    return result_list;
}

void GiDPostBinaryReader::ReadGaussPointRecord(const std::string& GpName)
{
    Reset();

    for(std::size_t i = 0; i < mGpRecord.size(); ++i)
    {
        bool check = (mGpRecord[i].Name.compare(GpName) == 0);
        check = (mGpRecord[i].NaturalCoordinates == "Given");

        if(!check)
            continue;

        /* iterate through gp coordinates */
        SetCurrentPosition(mGpRecord[i].StartCoordPos);

        std::vector<double> values;
        gid_value_t v;
        z_off_t CurPos;
        while(true)
        {
            Read(v);

            CurPos = GetCurrentPosition();
            values.push_back(static_cast<double>(v));

            if(CurPos >= mGpRecord[i].EndCoordPos)
            {
                break;
            }
        }
        values.pop_back(); // remove the last element. For some reason, GidPost library appends a dummy value at the end of the Gp coordinates.

        int dim = values.size() / mGpRecord[i].NumberOfGaussPoints;

        mGpRecord[i].GpCoordinates.clear();
        for (int j = 0; j < mGpRecord[i].NumberOfGaussPoints; ++j)
        {
            std::vector<double> p(dim);
            for (int k = 0; k < dim; ++k)
                p[k] = values[dim*j + k];
            mGpRecord[i].GpCoordinates.push_back(p);
        }
    }
}

void GiDPostBinaryReader::GetGaussPointRecordInfo(const std::string& GpName, int& Np, std::string& ElemType) const
{
    for(std::size_t i = 0; i < mGpRecord.size(); ++i)
    {
        if (mGpRecord[i].Name.compare(GpName) == 0)
        {
            Np = mGpRecord[i].NumberOfGaussPoints;
            ElemType = mGpRecord[i].ElemType;
            return;
        }
    }

    KRATOS_ERROR << "Gauss point record " << GpName << " is not found";
}

void GiDPostBinaryReader::GetGaussPointRecordInfo(const std::string& GpName, int& Np, std::string& ElemType, std::string& NaturalCoordinates) const
{
    for(std::size_t i = 0; i < mGpRecord.size(); ++i)
    {
        if (mGpRecord[i].Name.compare(GpName) == 0)
        {
            Np = mGpRecord[i].NumberOfGaussPoints;
            ElemType = mGpRecord[i].ElemType;
            NaturalCoordinates = mGpRecord[i].NaturalCoordinates;
            return;
        }
    }

    KRATOS_ERROR << "Gauss point record " << GpName << " is not found";
}

void GiDPostBinaryReader::GetGaussPointRecordCoordinates(const std::string& GpName, std::vector<std::vector<double> >& rCoordinates) const
{
    for(std::size_t i = 0; i < mGpRecord.size(); ++i)
    {
        if (mGpRecord[i].Name.compare(GpName) == 0)
        {
            rCoordinates = mGpRecord[i].GpCoordinates;
        }
    }
}

void GiDPostBinaryReader::GetGaussPointRecordCoordinates(const std::string& GpName, std::vector<array_1d<double, 3> >& rCoordinates) const
{
    for(std::size_t i = 0; i < mGpRecord.size(); ++i)
    {
        if (mGpRecord[i].Name.compare(GpName) == 0)
        {
            rCoordinates.resize(mGpRecord[i].GpCoordinates.size());
            for (std::size_t j = 0; j < mGpRecord[i].GpCoordinates.size(); ++j)
            {
                unsigned int dim = mGpRecord[i].GpCoordinates[j].size() < 3 ? mGpRecord[i].GpCoordinates[j].size() : 3;
                rCoordinates[j].clear();
                for (unsigned int k = 0; k < dim; ++k)
                    rCoordinates[j][k] = mGpRecord[i].GpCoordinates[j][k];
            }
        }
    }
}

void GiDPostBinaryReader::ReadGaussPointScalarValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<double> > >& rValues)
{
    Reset();

    step_list.clear();
    rValues.clear();

    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        std::string RefName = '"' + Name + '"';
        std::string RefGpName = '"' + GpName + '"';

        bool check = (mPosResult[i].Name.compare(RefName) == 0) || (mPosResult[i].Name.compare(Name) == 0);
        check = check && (mPosResult[i].Type.compare("Scalar") == 0);
        check = check && (mPosResult[i].Location.compare("OnGaussPoints") == 0);
        check = check && ((mPosResult[i].GpName.compare(RefGpName) == 0) || (mPosResult[i].GpName.compare(GpName) == 0));
        int np;
        bool found = false;
        for(std::size_t j = 0; j < mGpRecord.size(); ++j)
        {
            if(mGpRecord[j].Name == GpName || mGpRecord[j].Name == RefGpName)
            {
                found = true;
                np = mGpRecord[j].NumberOfGaussPoints;
            }
        }
        if(!found)
        {
            ERROR2("The GpRecord", RefGpName, "does not exist")
        }

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
                        std::vector<double> gp_values(np);
                        for(int j = 0; j < np; ++j)
                        {
                            Read(v);
                            gp_values[j] = v;
                        }
                        if(id > 0)
                        {
                            rValues[id].push_back(gp_values);
                            // REMARKS: for some reason, the last line contains the value pair: [id, v] = [-1, ##]. Probably to terminate the sequence of value. I don't know how to exclude this when I compute the end cursor position, so I make a naive comparison to make sure correct values are filtered out.
                        }
                        CurPos = GetCurrentPosition();
                    }
                    while(CurPos < mPosResult[i].EndPos);
                }
                else
                {
                    for(z_off_t i = 0; i < indexed; ++i)
                    {
                        Read(id);
                        std::vector<double> gp_values(np);
                        for(int j = 0; j < np; ++j)
                        {
                            Read(v);
                            gp_values[j] = v;
                        }
                        rValues[id].push_back(gp_values);
                    }
                }
            }
        }
    }
}

void GiDPostBinaryReader::ReadGaussPointVectorValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<std::vector<double> > > >& rValues)
{
    // TODO
}

void GiDPostBinaryReader::ReadGaussPointMatrixValues(const std::string& Name, const std::string& GpName, std::vector<double>& step_list, std::map<std::size_t, std::vector<std::vector<std::vector<double> > > >& rValues)
{
    Reset();

    step_list.clear();
    rValues.clear();

    for(std::size_t i = 0; i < mPosResult.size(); ++i)
    {
        std::string RefName = '"' + Name + '"';
        std::string RefGpName = '"' + GpName + '"';

        bool check = (mPosResult[i].Name.compare(RefName) == 0) || (mPosResult[i].Name.compare(Name) == 0);
        check = check && (mPosResult[i].Type.compare("Matrix") == 0);
        check = check && (mPosResult[i].Location.compare("OnGaussPoints") == 0);
        check = check && ((mPosResult[i].GpName.compare(RefGpName) == 0) || (mPosResult[i].GpName.compare(GpName) == 0));
        int np;
        bool found = false;
        for(std::size_t j = 0; j < mGpRecord.size(); ++j)
        {
            if(mGpRecord[j].Name == GpName || mGpRecord[j].Name == RefGpName)
            {
                found = true;
                np = mGpRecord[j].NumberOfGaussPoints;
            }
        }
        if(!found)
        {
            ERROR2("The GpRecord", RefGpName, "does not exist")
        }

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
                        std::vector<std::vector<double> > gp_values(np);
                        for(int j = 0; j < np; ++j)
                        {
                            gp_values[j].resize(6); // matrix result is alwayss symmetric 3x3 matrix stored as vector of size 6
                            for(int k = 0; k < 6; ++k)
                            {
                                Read(v);
                                gp_values[j][k] = v;
                            }
                        }
                        if(id > 0)
                        {
                            rValues[id].push_back(gp_values);
                            // REMARKS: for some reason, the last line contains the value pair: [id, v] = [-1, ##]. Probably to terminate the sequence of value. I don't know how to exclude this when I compute the end cursor position, so I make a naive comparison to make sure correct values are filtered out.
                        }
                        CurPos = GetCurrentPosition();
                    }
                    while(CurPos < mPosResult[i].EndPos);
                }
                else
                {
                    for(z_off_t i = 0; i < indexed; ++i)
                    {
                        Read(id);
                        std::vector<std::vector<double> > gp_values(np);
                        for(int j = 0; j < np; ++j)
                        {
                            gp_values[j].resize(6); // matrix result is alwayss symmetric 3x3 matrix stored as vector of size 6
                            for(int k = 0; k < 6; ++k)
                            {
                                Read(v);
                                gp_values[j][k] = v;
                            }
                        }
                        rValues[id].push_back(gp_values);
                    }
                }
            }
        }
    }
}

void GiDPostBinaryReader::ReadMesh(const std::string& Name, std::map<int, std::vector<double> >& rCoordinates)
{
    Reset();

    std::string RefName = '"' + Name + '"';
    std::vector<double> v(3);
    std::vector<gid_index_t> e;
    for(std::size_t i = 0; i < mPosMesh.size(); ++i)
    {
        bool check = (mPosMesh[i].Name.compare(RefName) == 0) || (mPosMesh[i].Name.compare(Name) == 0);

        if(!check)
            continue;

        bool coordinates_are_read = false;
        std::string line;
        std::vector<std::string> fields;
        z_off_t CurPos;

        if(mPosMesh[i].HasCoord())
        {
            /* iterate through all the coordinates */
            SetCurrentPosition(mPosMesh[i].StartCoordPos);

            do
            {
                ReadString(line);
//                        WATCH(line)

                fields.clear();
                boost::split(fields, line, boost::is_any_of(" "));

                if(fields[0].compare("Coordinates") == 0)
                {
                    z_off_t last_node_id = atoi(fields[1].c_str());
                    gid_index_t id;

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
        } // end if(mPosMesh[i].HasCoord())
    }
}

void GiDPostBinaryReader::ReadMesh(const std::string& Name, std::map<int, std::vector<int> >& rConnectivities)
{
    Reset();

    std::string RefName = '"' + Name + '"';
    std::vector<double> v(3);
    std::vector<gid_index_t> e;
    for(std::size_t i = 0; i < mPosMesh.size(); ++i)
    {
        bool check = (mPosMesh[i].Name.compare(RefName) == 0) || (mPosMesh[i].Name.compare(Name) == 0);

        if(!check)
            continue;

        bool connectivities_are_read = false;
        std::string line;
        std::vector<std::string> fields;
        z_off_t CurPos;

        if(mPosMesh[i].HasElem())
        {
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
                            if(e.size() != mPosMesh[i].Nnode)
                            {
                                e.resize(mPosMesh[i].Nnode);
                            }
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
        } // if(mPosMesh[i].HasElem())
    }
}

void GiDPostBinaryReader::ParseAndIndex2()
{
    auto time_begin = std::chrono::high_resolution_clock::now();

    char c;
    std::map<std::size_t, std::vector<z_off_t> > index;
    z_off_t temp_curr;
    bool found, next;
    int error;

    std::vector<std::string> keywords;
    keywords.push_back("MESH");
    keywords.push_back("Result");
    keywords.push_back("GaussPoints");

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

    //// DEBUGGING
    // for(std::size_t i = 0; i < keywords.size(); ++i)
    // {
    //     std::cout << "Indexing for keyword " << keywords[i] << ":";
    //     for(std::size_t j = 0; j < index[i].size(); ++j)
    //         std::cout << " " << index[i][j];
    //     std::cout << std::endl;
    // }
    //// END OF DEBUGGING

    // handle keyword "MESH"
    for(std::size_t i = 0; i < index[0].size(); ++i)
    {
        // set the cursor position
        SetCurrentPosition(index[0][i] + 5);

        // read the header
        MeshHeader mesh;
        ReadWordWithQuote(mesh.Name);

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

        mesh.StartCoordPos = mesh.EndCoordPos = mesh.StartElemPos = mesh.EndElemPos = 0;

        temp_curr = GetCurrentPosition();

        char first_str[80];
        error = ReadString(first_str, 80);
        if(error != 0) {next = false;} else {next = !strcmp(first_str, "Coordinates -1 Indexed");}
        if(next)
        {
            // search for the next End keyword
            found = false;
            std::string end_keyword = "End Coordinates";
            while(!CheckEof())
            {
                found = FindNext(0, c, end_keyword.c_str(), end_keyword.size());
                if(found)
                    break;
            }

            if(found)
            {
                mesh.StartCoordPos = temp_curr;
                mesh.EndCoordPos = GetCurrentPosition() - end_keyword.size()*sizeof(char);
            }
            else
            {
                ERROR("That should not happen", "No matching End Coordinates keyword")
            }

            // read through the coordinates
            std::string dummy;
            ReadWord(dummy);
        }
        else
            SetCurrentPosition(temp_curr);

        temp_curr = GetCurrentPosition();

        error = ReadString(first_str, 80);
        if(error != 0) {next = false;} else {next = !strcmp(first_str, "Elements -1 Indexed");}
        if(next)
        {
            // search for the next End keyword
            found = false;
            std::string end_keyword = "End Elements";
            while(!CheckEof())
            {
                found = FindNext(0, c, end_keyword.c_str(), end_keyword.size());
                if(found)
                    break;
            }
            if(found)
            {
                mesh.StartElemPos = temp_curr;
                mesh.EndElemPos = GetCurrentPosition() - end_keyword.size()*sizeof(char);
            }
            else
            {
                ERROR("That should not happen", "No matching End Elements keyword")
            }
        }
        else
            SetCurrentPosition(temp_curr);

        mPosMesh.push_back(mesh);
    } // end handle keyword "MESH"

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

        if(result.Location == std::string("OnGaussPoints"))
        {
            ReadWord(result.GpName);
        }

        result.StartPos = result.EndPos = 0;

        temp_curr = GetCurrentPosition();

        char first_str[80];
        error = ReadString(first_str, 80);
        if(error != 0) {next = false;} else {next = !strcmp(first_str, "Values -1 Indexed");}
        if(next)
        {
            // search for the next End keyword
            found = false;
            std::string end_keyword = "End Values";
            while(!CheckEof())
            {
                found = FindNext(0, c, end_keyword.c_str(), end_keyword.size());
                if(found)
                    break;
            }
            if(found)
            {
                result.StartPos = temp_curr;
                result.EndPos = GetCurrentPosition() - end_keyword.size()*sizeof(char);
            }
            else
            {
                ERROR("That should not happen", "No matching End Values keyword")
            }
        }
        else
            SetCurrentPosition(temp_curr);

        mPosResult.push_back(result);
    } // end handle keyword "Result"

    // handle keyword "GaussPoints"
    for(std::size_t i = 0; i < index[2].size(); ++i)
    {
        // set the cursor position
        SetCurrentPosition(index[2][i] + 12);

        // read the header
        GpRecord gp_rec;
        ReadWord(gp_rec.Name);

        // read the element type
        std::string dummy;
        ReadWord(dummy);

        if(dummy == "ElemType")
        {
            ReadWord(gp_rec.ElemType);

            // read the number of Gauss points
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            dummy.clear(); ReadWord(dummy);
            gp_rec.NumberOfGaussPoints = atoi(dummy.c_str());

            dummy.clear(); ReadWord(dummy); //WATCH(dummy)
            dummy.clear(); ReadWord(dummy); //WATCH(dummy)
            dummy.clear(); ReadWord(dummy); //WATCH(dummy)
            dummy.clear(); ReadWord(dummy); //WATCH(dummy)
            dummy.clear(); ReadWord(dummy); //WATCH(dummy)
            ReadWord(gp_rec.NaturalCoordinates);

            gp_rec.StartCoordPos = gp_rec.EndCoordPos = 0;

            temp_curr = GetCurrentPosition();

            next = true;
            if(next)
            {
                // search for the next End keyword
                found = false;
                std::string end_keyword = "End GaussPoints";
                while(!CheckEof())
                {
                    found = FindNext(0, c, end_keyword.c_str(), end_keyword.size());
                    if(found)
                        break;
                }
                if(found)
                {
                    gp_rec.StartCoordPos = temp_curr;
                    gp_rec.EndCoordPos = GetCurrentPosition() - end_keyword.size()*sizeof(char);
                }
                else
                {
                    ERROR("That should not happen", "No matching End GaussPoints keyword")
                }
            }
            else
                SetCurrentPosition(temp_curr);

            mGpRecord.push_back(gp_rec);
        }
    } // end handle keyword "GaussPoints"

    auto time_end = std::chrono::high_resolution_clock::now();
    auto elapsed = time_end - time_begin;
    long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    double milliseconds = microseconds / 1e3;
    std::cout << "Indexing " + BaseType::GetFileName() + " completed... " << milliseconds << "ms" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int GiDPostBinaryReader::Close()
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
int GiDPostBinaryReader::Open( const char * name )
{
    Close();
    m_file = gzopen(name, "rb");
    return m_file == NULL;
}

int GiDPostBinaryReader::Open( const std::string name )
{
    return Open(name.c_str());
}

/* initialize reading operation */
int GiDPostBinaryReader::Init()
{
    /* skip the first int */
    return gzseek(m_file, sizeof(int), SEEK_SET);
}

/* Reset the cursor position */
int GiDPostBinaryReader::Reset()
{
    return gzrewind(m_file);
}

/* read a string in GiD binary post file */
/* this version does not return error code */
//    char* ReadString(char* str)
//    {
//        int len;
//        gzread(m_file, &len, sizeof(int));
//        if(strlen(str) < len)
//            str = (char*)malloc(sizeof(char) * len);
//        gzread(m_file, str, len*sizeof(char));
//        return str;
//    }
void GiDPostBinaryReader::ReadString(char* str) const
{
    int len;
    gzread(m_file, &len, sizeof(int));
    WATCH(len)
    gzread(m_file, str, len*sizeof(char));
    str[len] = '\0';
}

int GiDPostBinaryReader::ReadString(char* str, int max_len) const
{
    int len;
    gzread(m_file, &len, sizeof(int));
    if((len > max_len) || (len < 0))
        return 1;
    gzread(m_file, str, len*sizeof(char));
    str[len] = '\0';
    return 0;
}

//char* GiDPostBinaryReader::ReadString(char* str)
//{
//    int len;
//    gzread(m_file, &len, sizeof(int));
//    if(strlen(str) < len)
//        str = (char*)malloc(sizeof(char) * len);
//    gzread(m_file, str, len*sizeof(char));
//    return str;
//}

int GiDPostBinaryReader::ReadString(std::string& str) const
{
    int len;
    int error_code = gzread(m_file, &len, sizeof(int)); // read the length of the string
    GZ_WATCH(error_code);
    if(len > 0)
    {
        if(str.size() != len)
            str.resize(len);
        error_code = gzread(m_file, (char*) str.c_str(), len*sizeof(char));
        GZ_WATCH(error_code);
        return 0;
    }
    else
    {
        str.resize(0);
        return 1;
    }
}

int GiDPostBinaryReader::CheckEof() const
{
    return gzeof(m_file);
}

z_off_t GiDPostBinaryReader::GetCurrentPosition() const
{
    return gztell(m_file);
}

z_off_t GiDPostBinaryReader::SetCurrentPosition(const z_off_t pos) const
{
    return gzseek(m_file, pos, SEEK_SET);
}

void GiDPostBinaryReader::ReadWord(std::string& word) const
{
    char c;
    Read(c);
    while((c != ' ') && (c != '\0'))
    {
        word.push_back(c);
        Read(c);
    }
}

void GiDPostBinaryReader::ReadWordWithQuote(std::string& word) const
{
    char c;
    Read(c);
    char c1 = c;
    int cntq = (c1 == '"') ? 1 : 0;
    while(((c != ' ') && (c != '\0'))
       || (c1 == '"' && cntq < 2))
    {
        word.push_back(c);
        Read(c);
        if (c == '"') ++cntq;
    }
}

/* shift n bytes */
z_off_t GiDPostBinaryReader::Shift(const z_off_t n) const
{
    return gzseek(m_file, n, SEEK_CUR);
}

bool GiDPostBinaryReader::isStringChar(const char ch) const {
    if (ch >= 'a' && ch <= 'z')
        return true;
    if (ch >= 'A' && ch <= 'Z')
        return true;
    return false;
}

bool GiDPostBinaryReader::isCharInString(const char c, const std::string& str) const
{
    for(int i = 0; i < str.size(); ++i)
    {
        if(c == str[i])
            return true;
    }
    return false;
}

bool GiDPostBinaryReader::CheckForString(const std::string& str) const
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

bool GiDPostBinaryReader::FindNext(int pos, char& c, const char* str, unsigned int len) const
{
    if(pos == len)
    {
        return true;
    }
    else
    {
        Read(c);
        if(str[pos] == ' ')
        {
            if(c == str[pos])
                return FindNext(pos+1, c, str, len);
            else
                return false;
        }
        else
        {
            if(c == tolower(str[pos]) || c == toupper(str[pos]))
                return FindNext(pos+1, c, str, len);
            else
                return false;
        }
    }
}

} // Namespace Kratos
