#include "MyMainWindow.h"
#include "OcctQtViewer.h"
#include"FileDialog.h"
#include <AIS_ViewCube.hxx>
#include <asi.hxx>
#include <occ.hxx>
#include <filesystem>
#include <common.hxx>

MyMainWindow::MyMainWindow() : myViewer(nullptr)
{
    setmenubar();
    setcheckboxes();
}

void MyMainWindow::setmenubar() {

    // menu bar with Quit item
    QMenuBar* aMenuBar = new QMenuBar();
    QMenu* aMenuWindow = aMenuBar->addMenu("&File");

    QAction* anActionLoad = new QAction(aMenuWindow);
    anActionLoad->setText("&Load");
    aMenuWindow->addAction(anActionLoad);
    connect(anActionLoad, &QAction::triggered, [this]() {loadPath(); });
    
    QMenu* computeMenu = aMenuWindow->addMenu("&Compute");

    QAction* anActionVMC = new QAction(computeMenu);
    anActionVMC->setText("VMC");
    computeMenu->addAction(anActionVMC);
    connect(anActionVMC, &QAction::triggered, [this](){compute();});
   
    QAction* anActionCNC = new QAction(computeMenu);
    anActionCNC->setText("CNC");
    computeMenu->addAction(anActionCNC);
    connect(anActionCNC, &QAction::triggered, [this](){compute();});

    QAction* anActionLaser = new QAction(computeMenu);
    anActionLaser->setText("Laser");
    computeMenu->addAction(anActionLaser);
    connect(anActionLaser, &QAction::triggered, [this]() {compute(); });

    QAction* anActionQuit = new QAction(aMenuWindow);
    anActionQuit->setText("Quit");
    aMenuWindow->addAction(anActionQuit);
    connect(anActionQuit, &QAction::triggered, [this](){close();});

    setMenuBar(aMenuBar);
}

void MyMainWindow::setcheckboxes() {

    //3D Viewer and some controls on top of it
    //Checkbox Group layout
    myViewer = new OcctQtViewer();
    QVBoxLayout* aLayout = new QVBoxLayout(myViewer);
    aLayout->setDirection(QBoxLayout::TopToBottom);
    aLayout->setAlignment(Qt::AlignTop);

    QButtonGroup* group = new QButtonGroup(myViewer);
    group->addButton(&wireframeCheckBox);
    group->addButton(&edgeCheckBox);
    group->addButton(&verticesCheckBox);
    group->setExclusive(false);// Mutual Exclusive

    //Checkbox
    //WireFrame checkbox
    wireframeCheckBox.setText("WireFrame");
    wireframeCheckBox.setStyleSheet("QCheckBox { color : white; }; QCheckBox::indicator { color:black; }");
    aLayout->addWidget(&wireframeCheckBox);
    wireframeCheckBox.move(30, 60);
    connect(&wireframeCheckBox, &QCheckBox::stateChanged, [&] {checkAndView(); });

    //edges checkbox
    edgeCheckBox.setText("Edges");
    edgeCheckBox.setStyleSheet("QCheckBox { color : white; }; QCheckBox::indicator { color:black; }");
    aLayout->addWidget(&edgeCheckBox);
    connect(&edgeCheckBox, &QCheckBox::stateChanged, [&] {checkAndView(); });

    //vertices checkbox
    verticesCheckBox.setText("Vertices");
    verticesCheckBox.setStyleSheet("QCheckBox { color : white; }; QCheckBox::indicator { color:black; }");
    aLayout->addWidget(&verticesCheckBox);
    connect(&verticesCheckBox, &QCheckBox::stateChanged, [&] {checkAndView(); });

    setCentralWidget(myViewer);
}

void MyMainWindow::firstLoadView()
{
    aShape.Nullify();
    myViewer->Context()->RemoveAll(true);

    //setting cube
    Handle(AIS_ViewCube) myViewCube = new AIS_ViewCube();
    myViewCube->SetFixedAnimationLoop(false);
    myViewCube->SetAutoStartAnimation(true);
    myViewCube->TransformPersistence()->SetOffset2d(Graphic3d_Vec2i(100, 150));
    myViewer->Context()->Display(myViewCube, 0, 0, false);
    aShape = new AIS_Shape(shape);
    myViewer->Context()->Display(aShape, AIS_Shaded, 0, true);
}

void MyMainWindow::onlyVerticesView()
{
    TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);

    // Extract vertices
    arryOFPoints = new TColgp_HArray1OfPnt(1, vertices.Extent());
    for (Standard_Integer i = 1; i <= vertices.Extent(); i++)
    {
        TopoDS_Vertex vertex = TopoDS::Vertex(vertices.FindKey(i));
        gp_Pnt currentGeometricPoint = BRep_Tool::Pnt(vertex);
        arryOFPoints->SetValue(i, currentGeometricPoint);
        aisShapeVertices = new AIS_Shape(vertex);
        aisShapeVertices->SetDisplayMode(AIS_DisplayMode::AIS_Shaded);
        myViewer->Context()->Display(aisShapeVertices, true);

        // After display, delete the aisShapeVertices object to free memory
        aisShapeVertices.Nullify();
    }
    
    // Now, clear the map after use (if it's a dynamic container)
    vertices.Clear();

    // Finally, delete arryOFPoints to avoid memory leak
    arryOFPoints.Nullify();
}

void MyMainWindow::onlyEdgeView(){
    //Shaded
    aShape = new AIS_Shape(shape);
    aShape->SetColor(defaultcol);
    aShape->SetDisplayMode(AIS_DisplayMode::AIS_Shaded);

    //WireFrame
    aShape = new AIS_Shape(shape);
    aShape->SetColor(Quantity_Color(0.0, 1.0, 1.0, Quantity_TOC_RGB));
    aShape->SetDisplayMode(AIS_DisplayMode::AIS_WireFrame);
    myViewer->Context()->Display(aShape, true);
}

void MyMainWindow::loadPath()
{
    aShape.Nullify();
    myViewer->Context()->RemoveAll(true);
    if(!CADFilePath.isEmpty())
    CADFilePath.clear();
    //setting cube
    Handle(AIS_ViewCube) myViewCube = new AIS_ViewCube();
    myViewCube->SetFixedAnimationLoop(false);
    myViewCube->SetAutoStartAnimation(true);
    myViewCube->TransformPersistence()->SetOffset2d(Graphic3d_Vec2i(100, 150));
    myViewer->Context()->Display(myViewCube, 0, 0, false);

    wireframeCheckBox.setCheckState(Qt::CheckState::Unchecked);
    verticesCheckBox.setCheckState(Qt::CheckState::Unchecked);
    edgeCheckBox.setCheckState(Qt::CheckState::Unchecked);
    

    try{
        FileDialog fileDialog;
        fileDialog.exec();

        auto path = fileDialog.selectedFiles();
        // std::cout << "filename:" << path.at(0).toStdString() << std::endl;
        QFileInfo fileInfo(path.at(0));
        if (fileInfo.isFile()) {
            IsloadedPathISFile = true;
            CADFilePath.insert(0, path.at(0));
            // QMessageBox::information(this, "Selected File", "You selected a file: " + path.at(0));
            if (!path.isEmpty())
            {

                // STEP file reader
                STEPControl_Reader reader;
                auto status = reader.ReadFile(path.at(0).toStdString().c_str());


                if (status == IFSelect_RetDone)
                {

                    // Load all the shapes in the STEP file
                    reader.TransferRoots();

                    // Retrieve the first shape from the STEP file
                    shape = reader.OneShape();


                    // You can now work with the shape, such as saving it as a BRep file:
                    BRepTools::Write(shape, "output.brep");
                    std::cout << shape.ShapeType() << std::endl;

                    // TopoDS_Shape aBox = BRepPrimAPI_MakeBox(100.0, 50.0, 90.0).Shape();
                    aShape = new AIS_Shape(shape);
                    aShape->Color(defaultcol);
                    myViewer->Context()->Display(aShape, AIS_Shaded, 0, true);

                    std::cout << "STEP file loaded successfully!" << std::endl;
                }
                else
                {
                    qDebug() << "Failed to read STEP file.";
                }
            }

        }
        else if (fileInfo.isDir()) {
            IsloadedPathISFile = false;
            QMessageBox::information(this, "Selected Directory", "You selected a directory: " + path.at(0));
        }
        fileDialog.close();
    }
    catch(...){
        qDebug() << "ERROR!!!";
    }

    myViewer->View()->Invalidate();
    myViewer->update();
}

void MyMainWindow::compute()
{
    try {
        if (!aShape.IsNull() || !shape.IsNull()) {
            QString filePath = QFileDialog::getSaveFileName(this, "Save Compute File", "D:\\CADFeatureRec\\test CAD models\\comput output", "Text Files (*.txt)");

            if (outFile.is_open())
                {
                    static int iteration = 0;

                    auto singleasObj = new AS_FeatureRecognition(CADFilePath.toUtf8().constData());

                    // Redirect cout to write to the file
                    std::streambuf* originalBuf = std::cout.rdbuf();  // Save the original buffer
                    if(CONSOLE)
                    std::cout.rdbuf(outFile.rdbuf());  // Redirect cout to the file

                    singleasObj->RunFeatureRecognition(CADFilePath.toUtf8().constData(), iteration++);

                    if (CONSOLE)
                    // Reset the buffer to print to the console again
                    std::cout.rdbuf(originalBuf);
                    
                    // Close the file stream
                    outFile.close();

                    std::cout << "Compute file saved successfully! at :" << filePath.toStdString() << std::endl;
                }
                else {
                qDebug() << "file Path is empty.";
            }
        }
        else {
            qDebug() << "Input Shape is empty.";
        }

    }
    catch (...) {
        qDebug() << "ERROR!!!";
    }
}

void MyMainWindow::checkAndView(){
    auto isWireframeModeChecked = wireframeCheckBox.checkState() == Qt::CheckState::Checked;
    auto isVerticesModeChecked = verticesCheckBox.checkState() == Qt::CheckState::Checked;
    auto isEdgeModeChecked = edgeCheckBox.checkState() == Qt::CheckState::Checked;

    try {
        if (!aShape.IsNull() || !shape.IsNull()) {
            if (isWireframeModeChecked && isEdgeModeChecked && isVerticesModeChecked) {
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
                onlyVerticesView();
                onlyEdgeView();
            }
            else if (isWireframeModeChecked && isEdgeModeChecked && !isVerticesModeChecked)
            {
                aisShapeVertices.Nullify();
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
                onlyEdgeView();
            }
            else if (isWireframeModeChecked && !isEdgeModeChecked && isVerticesModeChecked)
            {
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
                onlyVerticesView();
            }
            else if (!isWireframeModeChecked && isEdgeModeChecked && isVerticesModeChecked) {
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                onlyVerticesView();
                onlyEdgeView();

            }
            else if (isWireframeModeChecked && !isEdgeModeChecked && !isVerticesModeChecked)
            {
                aisShapeVertices.Nullify();
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
            }
            else if (!isWireframeModeChecked && isEdgeModeChecked && !isVerticesModeChecked)
            {
                aisShapeVertices.Nullify();
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                onlyEdgeView();
            }
            else if (!isWireframeModeChecked && !isEdgeModeChecked && isVerticesModeChecked)
            {
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
                onlyVerticesView();
            }
            else {
                aShape.Nullify();
                myViewer->Context()->RemoveAll(true);

                firstLoadView();
            }
        }
        else {
            if (isWireframeModeChecked || isEdgeModeChecked || isVerticesModeChecked) {
                wireframeCheckBox.setCheckState(Qt::CheckState::Unchecked);
                verticesCheckBox.setCheckState(Qt::CheckState::Unchecked);
                edgeCheckBox.setCheckState(Qt::CheckState::Unchecked);
                qDebug() << "Input Shape is empty.";
            }   
        }
    }
    catch (...) {
        qDebug() << "ERROR!!!";
    }
}