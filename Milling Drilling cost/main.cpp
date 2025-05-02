#include <iostream>
// project
#include <qt_main_window.hxx>
// 3rd-party: QT
#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QtWidgets/QApplication>

#define QT_UI 1
#define OCC 0

#if !QT_UI
#include <asi.hxx>
#include <common.hxx>
#include <occ.hxx>
#endif

int main(int argc, char* argv[])
{
#if QT_UI
    QApplication aQApp(argc, argv);
    QCoreApplication::setApplicationName("Cost-Estimation-Tool");
    QCoreApplication::setOrganizationName("Neilsoft");
    QCoreApplication::setApplicationVersion(QString("1.0"));
    MyMainWindow aMainWindow;
    aMainWindow.resize(aMainWindow.sizeHint());
    aMainWindow.show();
    return aQApp.exec();
#else
    bool status = true;
    OSD_Timer timer;
    std::vector<std::array<float, 3>> slotDimensions;
    std::vector<std::pair<float, float>> drilledHoleDimensions;
    std::vector<std::pair<float, float>> tappedHoleDimensions;
    std::vector<std::pair<float, float>> handleDimensions;
    std::array<float, 3> boxDimensions{ 0.0f };
    timer.Start();

    const char* path = //"D:\\Data\\Z_Feature_Recognition\\IGES\\Sample2.IGS";
        "D:\\Data\\Z_Feature_Recognition\\z_m2-5-x-5mm-a-f-x-st2-2x12mm-hexagon-self-tapping-male-female-without-undercut-1.snapshot.1\\M2.5 x 5mm A-F x ST2-2 x 12mm Hexagon Self Tapping Male-Female without Undercut Threaded Pillar.STEP";
    std::vector<std::string> filePaths;
    auto numCADFiles = Utility::EnumerateCADFiles(path, filePaths);

    if (-1 == numCADFiles)
    {
        // file path
        std::unique_ptr<AS_FeatureRecognition> asObj = std::make_unique<AS_FeatureRecognition>(path);
#if OCC
        auto occObj = asObj->GetCADPart();
        status &= occObj->ReadFile() ? false : true;
        if (status) status &= occObj->GenerateTopologyGeometryMaps();
        //if (status) status &= occObj->WriteFile(false) ? false : true;
        if (status) status &= occObj->Compute();
#else
        status = asObj->RunFeatureRecognition(filePaths.at(0),
            slotDimensions, drilledHoleDimensions, tappedHoleDimensions, handleDimensions, boxDimensions);
#endif
    }
    else if (numCADFiles > 0)
    {
        // folder path
        std::unique_ptr<AS_FeatureRecognition[]> asObjs = std::make_unique<AS_FeatureRecognition[]>(numCADFiles);
        for (int16_t i = 0; i < numCADFiles; ++i)
            status &= asObjs[i].RunFeatureRecognition(filePaths.at(i),
                slotDimensions, drilledHoleDimensions, tappedHoleDimensions, handleDimensions, boxDimensions);
    }
    else
        status &= false;

    std::cout << "\n\n";
    timer.Stop();
    timer.Show();
    return !status;
#endif
}

// TBD:
// use the mode of helix radius/pitch and sort the helix dimensions according to their axes
