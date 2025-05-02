#include <filesystem>
// project
#include <common.hxx>


int16_t Utility::EnumerateCADFiles( const char* path,
                                    std::vector<std::string>& filePaths )
{
    namespace fs = std::filesystem;

    std::string& folderPath = std::string(path);
    int16_t numCADFiles = 0;
    if (fs::is_directory(folderPath))
    {
        for (const auto& entry : fs::directory_iterator(folderPath))
        {
            if (entry.is_directory())
                continue;
            const std::string& cadFilePath = fs::absolute(entry.path()).string();
            std::string& cadFileExtn = Utility::FileExtn(cadFilePath.c_str());
            std::transform(cadFileExtn.begin(), cadFileExtn.end(), cadFileExtn.begin(), ::tolower);
            if (!(cadFileExtn == "stp" || cadFileExtn == "step" || cadFileExtn == "igs" || cadFileExtn == "iges"))
                continue;
            filePaths.emplace_back(cadFilePath);
            ++numCADFiles;
        }
    }
    else
    {
        filePaths.emplace_back(std::string(path));
        return -1;
    }

    return numCADFiles;
}


void Utility::SortMapByValue(   const std::map<int16_t, std::pair<float, float>>& map,
                                std::vector<std::pair<int16_t, float>>& vector  )
{
    vector.reserve(map.size());
    for (const auto& itr : map)
        vector.emplace_back(itr.first, itr.second.second);      // std::make_pair
    // comparator using lambda function
    std::sort(vector.begin(), vector.end(), \
        [](std::pair<int16_t, float>& p1, std::pair<int16_t, float>& p2) { return p1.second < p2.second; });
}


float Utility::ComputeMode( const std::vector<float>& vector,
                            int16_t& index, int16_t& frequency  )
{
    float mode = 0.0f;

    // map keeps a track of element frequency
    std::map<float, int16_t> map;
    for (const auto& vitr : vector)
    {
        auto mitr = map.find(vitr);
        if (map.end() == mitr)
            map.emplace(vitr, 1);
        else
            ++mitr->second;
    }
    for (auto itr = map.begin(); itr != map.end(); ++itr)
    {
        // select the element & index with highest frequency
        if (itr->second > frequency)
        {
            mode = itr->first;
            frequency = itr->second;
            index = std::distance(vector.begin(), std::find(vector.begin(), vector.end(), mode));
        }
    }

    return mode;
}
