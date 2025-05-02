#include <algorithm>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#define DEBUG 0
#define CONSOLE 1
#define RAD2DEG(x) ((180 / M_PI) * x)
#define DEG2RAD(x) ((M_PI / 180) * x)

class Utility
{
public:
    inline static std::string FileExtn(const char* filePath, bool lowerCase = true)
    {
        auto fileExtn = std::string(filePath).substr(std::string(filePath).rfind(".") + 1);
        if (lowerCase)
            std::transform(fileExtn.begin(), fileExtn.end(), fileExtn.begin(), ::tolower);
        return fileExtn;
    }

    inline static std::string FileName(const char* filePath)
    {
        return std::string(filePath).substr(std::string(filePath).rfind("\\") + 1);
    }

    inline static std::string FolderPath(const char* filePath)
    {
        return std::string(filePath).substr(0, std::string(filePath).rfind("\\") + 1);
    }

    inline static std::string OutputPath(const char* filePath, bool stepFormat)
    {
        auto folderPath = Utility::FolderPath(filePath);
        auto fileName = Utility::FileName(filePath);
        return (folderPath + fileName + (stepFormat ? ".step" : ".iges"));
    }

    static int16_t EnumerateCADFiles(   const char* dirPath,
                                        std::vector<std::string>& filePaths );

    static void SortMapByValue( const std::map<int, std::pair<float, float>>& map,
                                std::vector<std::pair<int, float>>& vector  );
};
