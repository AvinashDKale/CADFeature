#pragma once
#include "OcctQtViewer.h"
//std
#include <iostream>
#include <fstream>
#include <Standard_WarningsDisable.hxx>
#include <Standard_WarningsRestore.hxx>
#include <Standard_Version.hxx>

//Qt
#include <QFileDialog>
#include <QtCore/QtDebug>
#include <QtCore/QCommandLineParser>
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
#include <QVBoxLayout>
#include <QWidget>
#include <QCheckBox>
#include <QApplication>
#include <QSurfaceFormat>
#include <QAction>
#include <QLabel>
#include <QMainWindow>
#include <QMenuBar>
#include <QMessageBox>
#include <QVBoxLayout>
#include<QTextEdit>
#include <QPushButton>
#include <QSlider>
#include <QFrame>
#include <QButtonGroup>
#include <QPushButton>
#include <QTreeWidget>
#include <QStringList>
#include <QModelIndex>
#include <QDir>
#include <QDebug>
#include<QListView>



//opengl
#include <OpenGl_GraphicDriver.hxx>
#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>
#include <gp_Pnt.hxx> 

//opencascade
#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <STEPCAFControl_Controller.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx> 
#include <TopoDS_Edge.hxx> 
#include <BRep_Builder.hxx>
#include <V3d_Viewer.hxx>
#include <AIS_Shape.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <V3d_View.hxx>
#include <AIS_InteractiveContext.hxx>


//! Main application window.
class MyMainWindow : public QMainWindow
{

    QListView* m_listView;
    QTreeView* m_treeView;
    QPushButton* m_btnOpen;
    QStringList m_selectedFiles;


    OcctQtViewer* myViewer;
    TopoDS_Shape shape;
    Handle(AIS_Shape) aShape;
    QFrame frame;
    QCheckBox testcheckBox{&frame};
    QCheckBox wireframeCheckBox{ &frame };
    QCheckBox edgeCheckBox{ &frame };
    QCheckBox verticesCheckBox{ &frame };

    QString CADFilePath="";
    bool IsloadedPathISFile;

    //vertices,shaded view
    Handle(AIS_Shape) aisShape1;
    Handle(AIS_Shape) aisShapeVertices;

    //colour
    Quantity_Color defaultcol; 
    Quantity_Color defaultColWire;

    //to extract vertices
    TopTools_IndexedMapOfShape vertices;
    Handle_TColgp_HArray1OfPnt arryOFPoints;

    //AS_FeatureRecognition singleasObj;

    void setmenubar();
    void setcheckboxes();

    void firstLoadView();
    void onlyVerticesView();
    void onlyEdgeView();

private slots:
    void checkAndView();
    void loadPath();
    void compute();
    
public:
    MyMainWindow();
};

