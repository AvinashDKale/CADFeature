#include "MyMainWindow.h"
#include "OcctQtViewer.h"
#include <AIS_ViewCube.hxx>


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
    anActionLoad->setText("Load");
    aMenuWindow->addAction(anActionLoad);
    connect(anActionLoad, &QAction::triggered, [this]() { loadFile(); });

    QAction* anActionCompute = new QAction(aMenuWindow);
    anActionCompute->setText("Compute");
    aMenuWindow->addAction(anActionCompute);
    connect(anActionCompute, &QAction::triggered, [this]() {compute(); });

    QAction* anActionQuit = new QAction(aMenuWindow);
    anActionQuit->setText("Quit");
    aMenuWindow->addAction(anActionQuit);
    connect(anActionQuit, &QAction::triggered, [this]()
        {
            close();
        });

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
    }
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

void MyMainWindow::loadFile()
{
    if (!aShape.IsNull() || shape.IsNull())
    {
        myViewer->Context()->Remove(aShape, true);
        myViewer->Context()->Remove(aisShapeVertices, true);
        shape.Nullify(); aShape.Nullify();
        aisShapeVertices.Nullify();
    }

    QString filePath = QFileDialog::getOpenFileName(this, "Open STEP File", "D:\\CADFeatureRec\\test CAD models", "STEP Files (*.step *.stp)");
    if (!filePath.isEmpty())
    {

        // STEP file reader
        STEPControl_Reader reader;
        auto status = reader.ReadFile(filePath.toStdString().c_str());


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

    myViewer->View()->Invalidate();
    myViewer->update();
}

void MyMainWindow::compute()
{
    try {
        if (!aShape.IsNull() || !shape.IsNull()) {
            QString filePath = QFileDialog::getSaveFileName(this, "Save Compute File", "", "Text Files (*.txt)");

            if (!filePath.isEmpty())
            {

                // Open the chosen file for writing text
                QFile file(filePath);
                if (file.open(QIODevice::WriteOnly | QIODevice::Text))
                {
                    QTextStream out(&file);

                    // Write the shape type to the file
                    out << "Shape Type: " << shape.ShapeType() << "\n";

                    // Optionally write additional information to the file
                    out << "STEP file loaded successfully!" << "\n";

                    // Close the file after writing
                    file.close();
                    std::cout << "Compute file saved successfully! at :" << filePath.toStdString() << std::endl;
                }
                else
                {
                    qDebug() << "Failed to open file for writing.";
                }
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
                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
                onlyVerticesView();
                onlyEdgeView();
            }
            else if (isWireframeModeChecked && isEdgeModeChecked && !isVerticesModeChecked)
            {
                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
                onlyEdgeView();
            }
            else if (isWireframeModeChecked && !isEdgeModeChecked && isVerticesModeChecked)
            {
                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
                onlyVerticesView();
            }
            else if (!isWireframeModeChecked && isEdgeModeChecked && isVerticesModeChecked) {
                firstLoadView();
                onlyVerticesView();
                onlyEdgeView();

            }
            else if (isWireframeModeChecked && !isEdgeModeChecked && !isVerticesModeChecked)
            {
                firstLoadView();
                myViewer->Context()->Display(aShape, AIS_WireFrame, 0, true);
            }
            else if (!isWireframeModeChecked && isEdgeModeChecked && !isVerticesModeChecked)
            {
                firstLoadView();
                onlyEdgeView();
            }
            else if (!isWireframeModeChecked && !isEdgeModeChecked && isVerticesModeChecked)
            {
                firstLoadView();
                onlyVerticesView();
            }
            else {
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


