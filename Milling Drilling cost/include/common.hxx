#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#define PRODUCTION 0
#define RAD2DEG(x) ((180 / M_PI) * x)
#define DEG2RAD(x) ((M_PI / 180) * x)

class Utility
{
public:
    inline static bool CONSOLE = false;
    inline static int16_t DEBUG_LEVEL = 1;
    inline static void ConstrainDebugLevel()
    {
        if (PRODUCTION) DEBUG_LEVEL = 0;
    }

    inline static std::string FileExtn(const char* filePath, bool lowerCase = true)
    {
        auto fileExtn = std::string(filePath).substr(std::string(filePath).rfind(".") + 1);
        if (lowerCase)
            std::transform(fileExtn.begin(), fileExtn.end(), fileExtn.begin(), ::tolower);
        return fileExtn;
    }

    inline static std::string FileName(const char* filePath)
    {
        // UNIX-path: forward slash, WINDOWS-path: backward slash
        char str[2]{ '\0' };
        if (std::string::npos != std::string(filePath).find('/'))
            str[0] = '/';
        else
            str[0] = '\\';
        return std::string(filePath).substr(std::string(filePath).rfind(str[0]) + 1);
    }

    inline static std::string FolderPath(const char* filePath)
    {
        // UNIX-path: forward slash, WINDOWS-path: backward slash
        char str[2]{ '\0' };
        if (std::string::npos != std::string(filePath).find('/'))
            str[0] = '/';
        else
            str[0] = '\\';
        return std::string(filePath).substr(0, std::string(filePath).rfind(str[0]) + 1);
    }

    inline static std::string OutputPath(const char* filePath, int16_t format)
    {
        auto folderPath = Utility::FolderPath(filePath);
        auto fileName = Utility::FileName(filePath);
        auto fileExtn = !format ? ".step" : (1 == format ? ".iges" : ".brep");
        return (folderPath + fileName + fileExtn);
    }

    template<typename T>
    static inline int16_t NonZeroContainerCount(T container)
    {
        auto nonZeroCount = std::count_if(container.begin(), container.end(), \
            [](auto elem) -> bool { return (elem > 0); });
        return nonZeroCount;
    }

    static int16_t EnumerateCADFiles(   const char* path,
                                        std::vector<std::string>& filePaths );

    static void SortMapByValue( const std::map<int16_t, std::pair<float, float>>& map,
                                std::vector<std::pair<int16_t, float>>& vector  );

    static float ComputeMode(   const std::vector<float>& vector,
                                int16_t& index, int16_t& frequency  );
};
