#include <array>
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>



// 3rd-party: QT
#include <QAction>
#include <QApplication>
#include <QButtonGroup>
#include <QCheckBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFrame>
#include <QGroupBox>
#include <QComboBox>
#include <QLineEdit>
#include <QLabel>
#include <QMainWindow>
#include <QMenuBar>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPushButton>
#include <QRadioButton>
#include <QSettings>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QWidget>
#include <QWidgetAction>
#include <QtCore/QCommandLineParser>
#include <QtCore/QtDebug>
#include <QtCore/QDir>
#include <QtCore/QSettings>
#include <QtCore/QTimer>
#include <QtCore/QTranslator>
#include <QtCore/QVersionNumber>
#include <QtGui/QOffscreenSurface>
#include <QtGui/QOpenGLContext>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QPixmap>
#include <QtWidgets/QApplication>



class ReadMachiningParameters {
private:
    int material{ -1 };
    float rawMaterialRate{ -1 };
    float cutterDiameter{ -1 };
    float roughMachiningDepth{ -1 };
    float finishingMachiningDepth{ -1 };
    float millingRate{ -1 };
    float turningRate{ -1 };
    float drillingRate{ -1 };
    float tappingRate{ -1 };
    const char* filename = "D:\\CADFeatureRec\\code\\proj\\Demo4\\maching_data.txt";
public:


    std::map<int, std::string> materialLegend = {
        {0, "Steel tough"},
        {1, "Mild Steel"},
        {2, "Cast Iron"},
        {3, "Alloy steel"},
        {4, "Stainless steel"},
        {5, "Bronze"},
        {6, "Aluminium"}
    };

    std::map<int, int> materialRateLegend = {
    {0, 75},
    {1, 65},
    {2, 40},
    {3, 70},
    {4, 210},
    {5, 450},
    {6, 290}
    };

    std::map<int, int> CutterDiameterLegend = {
    {0, 6},
    {1, 10},
    {2, 20},
    {3, 25},
    {4, 32},
    {5, 40},
    {6, 50},
    {7, 63},
    {8, 80},
    {9, 100},
    {10, 120},
    {11, 160},
    {12, 220},
    {13, 250},
    {14, 300}
    };

    ReadMachiningParameters() {


        std::ifstream file(filename);

        if (!file) {
            std::cerr << "Error opening MachingParameters file!" << std::endl;

        }

        char line[256]; // Buffer for reading lines
        bool parsingLegends = false;

        while (file.getline(line, sizeof(line))) {
            // Skip empty lines
            if (std::strlen(line) == 0) continue;

            // Detect "legends" section and stop parsing
            if (std::strstr(line, "legends")) {
                parsingLegends = true;
                continue;
            }
            if (parsingLegends) continue; // Ignore legends section

            // Tokenize line using `strtok`
            const char* key = std::strtok(line, ":");
            char* valueStr = std::strtok(nullptr, ":");

            if (!key || !valueStr) continue; // Skip invalid lines

            float value = std::atof(valueStr); // Convert value to double

            // Map the extracted key-value pair to the class
            if (std::strlen("Material") == std::strlen(key) && std::strstr(key, "Material")) {

                material = static_cast<int>(std::atof(valueStr));
            }

            else if (std::strlen("Raw Material Rate") == std::strlen(key) && std::strstr(key, "Raw Material Rate")) rawMaterialRate = value;

            else if (std::strlen("Cutter Diameter") == std::strlen(key) && std::strstr(key, "Cutter Diameter"))
            {
                cutterDiameter = value;
                //// Iterate through the map to find the key for the value
                //for (const auto& pair : CutterDiameterLegend) {
                //    if (pair.second == cutterDiameter) {
                //        cutterDiameter = pair.first;
                //        break;
                //    }
                //}
            }


            else if (std::strlen("Rough Machining Depth") == std::strlen(key) && std::strstr(key, "Rough Machining Depth")) roughMachiningDepth = value;
            else if (std::strlen("Finishing Machining Depth") == std::strlen(key) && std::strstr(key, "Finishing Machining Depth")) finishingMachiningDepth = value;
            else if (std::strlen("Milling Rate") == std::strlen(key) && std::strstr(key, "Milling Rate")) millingRate = value;
            else if (std::strlen("Turning Rate") == std::strlen(key) && std::strstr(key, "Turning Rate")) turningRate = value;
            else if (std::strlen("Drilling Rate") == std::strlen(key) && std::strstr(key, "Drilling Rate")) drillingRate = value;
            else if (std::strlen("Tapping Rate") == std::strlen(key) && std::strstr(key, "Tapping Rate")) tappingRate = value;
        }

        file.close();
    }

    // Getters
    int getMaterial() const { return material; }
    int getRawMaterialRate() const { return rawMaterialRate; }
    int getCutterDiameterLegendKey() {
        int key = -1; // If no key is found, this will remain -1

        // Iterate through the map to find the key for the value
        for (const auto& pair : CutterDiameterLegend) {
            if (pair.second == cutterDiameter) {
                key = pair.first;
                break;
            }
        }
        return key;
    }
    float getCutterDiameter() const { return cutterDiameter; }
    float getRoughMachiningDepth() const { return roughMachiningDepth; }
    float getFinishingMachiningDepth() const { return finishingMachiningDepth; }
    float getMillingRate() const { return millingRate; }
    float getTurningRate() const { return turningRate; }
    float getDrillingRate() const { return drillingRate; }
    float getTappingRate() const { return tappingRate; }

    // Setters
    void setMaterial(int value) { material =/*auto v = materialLegend.at(value); */ value; }
    void setRawMaterialRate(int value) { rawMaterialRate = value; }
    void setCutterDiameter(int value) { cutterDiameter = value; }
    void setRoughMachiningDepth(float value) { roughMachiningDepth = value; }
    void setFinishingMachiningDepth(float value) { finishingMachiningDepth = value; }
    void setMillingRate(float value) { millingRate = value; }
    void setTurningRate(float value) { turningRate = value; }
    void setDrillingRate(float value) { drillingRate = value; }
    void setTappingRate(float value) { tappingRate = value; }


    // save updated values to a file
    void saveMachingParameters() {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Error opening file for writing!" << std::endl;
            return;
        }

        file << "Material: " << material << "\n";
        file << "Raw Material Rate: " << rawMaterialRate << "\n";
        file << "Cutter Diameter: " << cutterDiameter << "\n";
        file << "Rough Machining Depth: " << roughMachiningDepth << "\n";
        file << "Finishing Machining Depth: " << finishingMachiningDepth << "\n";
        file << "Milling Rate: " << millingRate << "\n";
        file << "Turning Rate: " << turningRate << "\n";
        file << "Drilling Rate: " << drillingRate << "\n";
        file << "Tapping Rate: " << tappingRate << "\n";

        file << "\n";
        file << "legends: " << "\n";
        file << "Material : int" << "\n";
        file << "Steel tough : 0" << "\n";
        file << "Mild Steel : 1" << "\n";
        file << "Cast Iron : 2" << "\n";
        file << "Alloy steel : 3" << "\n";
        file << "Stainless steel : 4" << "\n";
        file << "Bronze : 5" << "\n";
        file << "Aluminium : 6" << "\n";


        /* file << "Cutter Diameter : int" << "\n";
         file << "6 : 0 " << "\n";
         file << "10 : 1" << "\n";
         file << "20 : 2" << "\n";
         file << "25 : 3" << "\n";
         file << "32 : 4" << "\n";
         file << "40 : 5" << "\n";
         file << "50 : 6" << "\n";
         file << "63 : 7" << "\n";
         file << "80 : 8" << "\n";
         file << "100 : 9" << "\n";
         file << "120 : 10" << "\n";
         file << "160 : 11" << "\n";
         file << "220 : 12" << "\n";
         file << "250 : 13" << "\n";
         file << "300 : 14" << "\n";*/

        file.close();
        std::cout << "Updated Machining Parameters saved successfully!\n";
    }

};

class CostEstimation : public QWidget
{
    Q_OBJECT
public:

    CostEstimation() = default;
    CostEstimation(
        //std::shared_ptr<OcctQtViewer> myQtViewer ,
        int16_t rawMaterial,
        float roughCuttingDepth,
        float finishCuttingDepth,
        int16_t cutterDiameter,
        int16_t rawMaterialRate,
        int16_t millingRate = 150,
        int16_t turningRate = 150,
        int16_t drillingRate = 20,
        int16_t tappingRate = 40);
    CostEstimation( // cad dimensions : asi/occ output
       // std::shared_ptr<OcctQtViewer> myQtViewer,
        const std::vector<std::array<float, 3>>& slotDimensions,
        const std::vector<std::pair<float, float>>& holeDimensions,
        const std::vector<std::pair<float, float>>& handleDimensions,
        const std::array<float, 3>& boxDimensions,
        // cost estimation : user input
        int16_t rawMaterial,
        float roughCuttingDepth,
        float finishCuttingDepth,
        int16_t cutterDiameter,
        int16_t rawMaterialRate,
        int16_t millingRate = 150,
        int16_t turningRate = 150,
        int16_t drillingRate = 20,
        int16_t tappingRate = 40);
    int noDrillTooShallow = 0;


    // non-virtual
    ~CostEstimation();

    float RawMaterialCost(std::array<float, 3>* boxDimensions);
    float MillingCost(std::vector<std::array<float, 3>>* slotDimensions);
    float DrillingTappingCost(std::vector<std::pair<float, float>>* holeDimensions);
    float TurningCost();

    int RunCostEstimation(std::vector<std::array<float, 3>>* slotDimensions = nullptr,
        std::vector<std::pair<float, float>>* holeDimensions = nullptr,
        std::vector<std::pair<float, float>>* handleDimensions = nullptr,
        std::array<float, 3>* boxDimensions = nullptr);

    std::map<int, std::string> materialLegend = {
    {0, "Steel tough"},
    {1, "Mild Steel"},
    {2, "Cast Iron"},
    {3, "Alloy steel"},
    {4, "Stainless steel"},
    {5, "Bronze"},
    {6, "Aluminium"}
    };


private:
    // CAD dimensions
    std::shared_ptr<std::vector<std::array<float, 3>>> m_Slots{ nullptr };
    std::shared_ptr<std::vector<std::pair<float, float>>> m_Holes{ nullptr }, m_Handles{ nullptr };
    float m_BoxVolume{ 0 };

    // Manufacturing cost
    float m_RoughMillingCuttingDepth{ 0.0f }, m_FinishMillingCuttingDepth{ 0.0f };
    float m_MaterialDensity{ 0.0f }, m_MaterialRate{ 0.0f };
    int16_t m_RawMaterial;
    int16_t m_MillingRate{ 0 }, m_TurningRate{ 0 }, m_DrillingRate{ 0 }, m_TappingRate{ 0 };
    int16_t m_RoughCuttingFeed{ 0 }, m_FinishCuttingFeed{ 0 }, m_CutterDiameter{ 0 };

public:
    enum class RawMaterial { TOUGH_STEEL, MILD_STEEL, CAST_IRON, ALLOY_STEEL, STAINLESS_STEEL, BRONZE, ALUMINIUM };
    enum class CutterDiameter
    {
        SIX, TEN, TWENTY, TWENTY_FIVE, THIRTY_TWO, FORTY, FIFTY, SIXTY_THREE, EIGHTY, ONE_HUNDRED,
        ONE_HUNDRED_TWENTY, ONE_HUNDRED_SIXTY, TWO_HUNDRED_TWENTY, TWO_HUNDRED_FIFTY, THREE_HUNDRED
    };

    
private:
    QWidget* drillCostForm = nullptr;
    void closeEvent(QCloseEvent* event);
    bool useSpecialTool = false;
};
