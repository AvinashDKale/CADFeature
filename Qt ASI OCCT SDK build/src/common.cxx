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


void Utility::SortMapByValue(   const std::map<int, std::pair<float, float>>& map,
                                std::vector<std::pair<int, float>>& vector  )
{
    vector.reserve(map.size());
    for (const auto& itr : map)
        vector.emplace_back(std::make_pair(itr.first, itr.second.second));
    // comparator using lambda function
    std::sort(vector.begin(), vector.end(), \
        [](std::pair<int, float>& p1, std::pair<int, float>& p2) { return p1.second < p2.second; });
}
