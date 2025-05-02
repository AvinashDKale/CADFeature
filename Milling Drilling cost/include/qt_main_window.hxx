#include <fstream>
#include <iostream>
// project
#include <qt.hxx>
#include <cost.hxx>
// 3rd-party: OCC
#include <AIS_Shape.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <Standard_WarningsDisable.hxx>
#include <Standard_WarningsRestore.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx> 
#include <TopoDS_Vertex.hxx> 
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <V3d_Viewer.hxx>
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


class MyMainWindow : public QMainWindow
{
private:
    std::shared_ptr<OcctQtViewer> myQtViewer{ nullptr };
    std::shared_ptr<TopoDS_Shape> topoShape{ nullptr };
    Handle(AIS_Shape) aisShape;

    QFrame frame;
    QCheckBox testcheckBox{ &frame };
    QCheckBox wireframeCheckBox{ &frame };
    QCheckBox edgeCheckBox{ &frame };
    QCheckBox verticesCheckBox{ &frame };


    QWidget* costEstimationForm = nullptr;
    QComboBox* rawMaterialComboBox;
    QComboBox* cutterDiameterComboBox;
    QLineEdit* roughMachiningDepthLineEdit;
    QLineEdit* rawMaterialRateLineEdit;
    QLineEdit* finishingMachiningDepthLineEdit;
    QLineEdit* millingRateLineEdit;
    QLineEdit* turningRateLineEdit;
    QLineEdit* drillingRateLineEdit;
    QLineEdit* tappingRateLineEdit;
    QLabel* resultLabel;
    void onEstimateButtonClicked();
    ReadMachiningParameters readerMachiningParameters;
    void onShowCostEstimationForm();
    void onSaveButtonClicked();

    // vertices, shaded view
    Handle(AIS_Shape) aisShapeVertices;
    // extract vertices
    std::shared_ptr<TopTools_IndexedMapOfShape> vertices{ nullptr };
    Handle_TColgp_HArray1OfPnt pointArray;

    QButtonGroup myButtonGroup;

    // color
    Quantity_Color defaultColor;
    Quantity_Color defaultColWire;

    void setMenubar();
    void setCheckboxes();

    void firstLoadView();
    void onlyVerticesView();
    void onlyEdgeView();

private slots:
    void checkAndView();
    void loadFile();
    void compute(const std::shared_ptr<TopoDS_Shape>& topoShape);
    void reset(bool flag = true);

public:
    MyMainWindow();
    void closeEvent(QCloseEvent* event);
};
