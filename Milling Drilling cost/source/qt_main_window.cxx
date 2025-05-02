// project
#include <asi.hxx>
#include <common.hxx>
#include <occ.hxx>
#include <qt_main_window.hxx>
// 3rd-party: OCC
#include <AIS_ViewCube.hxx>


MyMainWindow::MyMainWindow() : myQtViewer(nullptr)
{
    this->topoShape = std::make_shared<TopoDS_Shape>();
    this->vertices = std::make_shared<TopTools_IndexedMapOfShape>();
    this->setMenubar();
    this->setCheckboxes();
}


void MyMainWindow::setMenubar()
{
    QSettings settings("YourCompany", "YourApp");
    int savedId = settings.value("printOption", 2).toInt();  // default to 2 = File

    QMenuBar* aMenuBar = new QMenuBar();
    QMenu* aMenuWindow1 = aMenuBar->addMenu("&File");

    QAction* anActionLoad = new QAction("Load", aMenuWindow1);
    aMenuWindow1->addAction(anActionLoad);
    connect(anActionLoad, &QAction::triggered, [this]() {
        this->loadFile();
        });

    // Create radio buttons
    QGroupBox* groupBox = new QGroupBox("Print");
    QVBoxLayout* groupLayout = new QVBoxLayout(groupBox);
    QRadioButton* button1 = new QRadioButton("Console", groupBox);
    QRadioButton* button2 = new QRadioButton("File", groupBox);
    groupLayout->addWidget(button1);
    groupLayout->addWidget(button2);
    groupBox->setLayout(groupLayout);

    // Add buttons to group
    this->myButtonGroup.addButton(button1, 1);
    this->myButtonGroup.addButton(button2, 2);

    // Restore saved selection
    QAbstractButton* savedButton = this->myButtonGroup.button(savedId);
    if (savedButton) {
        savedButton->setChecked(true);
        Utility::CONSOLE = (savedId == 1); // set Utility::CONSOLE accordingly
    }

    // Output menu with group box as widget action
    QMenu* outputMenu = aMenuWindow1->addMenu("Output");
    QWidgetAction* groupAction = new QWidgetAction(outputMenu);
    groupAction->setDefaultWidget(groupBox);
    outputMenu->addAction(groupAction);

    // Handle button selection change
    connect(button1, &QRadioButton::clicked, this, [=]() {
        if (button1->isChecked()) {
            Utility::CONSOLE = true;
            QSettings settings("YourCompany", "YourApp");
            settings.setValue("printOption", 1);
        }
        });

    connect(button2, &QRadioButton::clicked, this, [=]() {
        if (button2->isChecked()) {
            Utility::CONSOLE = false;
            QSettings settings("YourCompany", "YourApp");
            settings.setValue("printOption", 2);
        }
        });


    Utility::ConstrainDebugLevel();

    QAction* anActionCompute = new QAction(aMenuWindow1);
    anActionCompute->setText("Compute");
    aMenuWindow1->addAction(anActionCompute);
    connect(anActionCompute, &QAction::triggered, [this]() { this->compute(this->topoShape); });

    QAction* anActionReset = new QAction(aMenuWindow1);
    anActionReset->setText("Reset");
    aMenuWindow1->addAction(anActionReset);
    connect(anActionReset, &QAction::triggered, [this]() { this->reset(); });

    QAction* anActionQuit = new QAction(aMenuWindow1);
    anActionQuit->setText("Quit");
    aMenuWindow1->addAction(anActionQuit);
    connect(anActionQuit, &QAction::triggered, [this]() { close(); });

    QMenu* aMenuWindow2 = aMenuBar->addMenu("&About");
    QAction* anActionAbout = new QAction(aMenuWindow2);
    anActionAbout->setText("Version");
    aMenuWindow2->addAction(anActionAbout);
    connect(anActionAbout, &QAction::triggered, [this]()
    {
        QMessageBox::information(0, "About Application", QString()
            + "Manufacturing-Cost-Estimation-CAD-Feature-Recognition\n\n"
            + "Organization: Neil Automation" "\n"
            + "Product ver: 1.0" "\n"
            + "ASi ver: 2024.1" "\n"
            + "OCCT ver: " OCC_VERSION_STRING_EXT "\n"
            + "Qt ver: " QT_VERSION_STR "\n\n"
            + "OpenGL info:\n"
            + this->myQtViewer->getGlInfo());
    });

    this->setMenuBar(aMenuBar);
}


void MyMainWindow::setCheckboxes()
{
    //3D Viewer and some controls on top of it
    //Checkbox Group layout
    this->myQtViewer = std::make_shared<OcctQtViewer>();
    QVBoxLayout* aLayout = new QVBoxLayout(this->myQtViewer.get());
    aLayout->setDirection(QBoxLayout::TopToBottom);
    aLayout->setAlignment(Qt::AlignTop);

    QButtonGroup* group = new QButtonGroup(this->myQtViewer.get());
    group->addButton(&this->wireframeCheckBox);
    group->addButton(&this->edgeCheckBox);
    group->addButton(&this->verticesCheckBox);
    group->setExclusive(false);// Mutual Exclusive

    //Checkbox
    //WireFrame checkbox
    this->wireframeCheckBox.setText("WireFrame");
    this->wireframeCheckBox.setStyleSheet("QCheckBox { color : white; }; QCheckBox::indicator { color : black; }");
    aLayout->addWidget(&this->wireframeCheckBox);
    this->wireframeCheckBox.move(30, 60);
    connect(&this->wireframeCheckBox, &QCheckBox::stateChanged, [&] { this->checkAndView(); });

    //edges checkbox
    this->edgeCheckBox.setText("Edges");
    this->edgeCheckBox.setStyleSheet("QCheckBox { color : white; }; QCheckBox::indicator { color:black; }");
    aLayout->addWidget(&this->edgeCheckBox);
    connect(&this->edgeCheckBox, &QCheckBox::stateChanged, [&] { this->checkAndView(); });

    //vertices checkbox
    this->verticesCheckBox.setText("Vertices");
    this->verticesCheckBox.setStyleSheet("QCheckBox { color : white; }; QCheckBox::indicator { color:black; }");
    aLayout->addWidget(&this->verticesCheckBox);
    connect(&this->verticesCheckBox, &QCheckBox::stateChanged, [&] { this->checkAndView(); });

    setCentralWidget(this->myQtViewer.get());
}


void MyMainWindow::firstLoadView()
{
    this->aisShape.Nullify();
    this->myQtViewer->Context()->RemoveAll(true);

    //setting cube
    Handle(AIS_ViewCube) aisViewCube = new AIS_ViewCube();
    aisViewCube->SetFixedAnimationLoop(false);
    aisViewCube->SetAutoStartAnimation(true);
    aisViewCube->TransformPersistence()->SetOffset2d(Graphic3d_Vec2i(100, 150));
    myQtViewer->Context()->Display(aisViewCube, 0, 0, false);
    this->aisShape = new AIS_Shape(*this->topoShape);
    this->myQtViewer->Context()->Display(this->aisShape, AIS_Shaded, 0, true);
}


void MyMainWindow::onlyVerticesView()
{
    TopExp::MapShapes(*this->topoShape, TopAbs_VERTEX, *this->vertices);

    // Extract vertices
    this->pointArray = new TColgp_HArray1OfPnt(1, this->vertices->Extent());
    for (Standard_Integer i = 1; i <= vertices->Extent(); i++)
    {
        TopoDS_Vertex vertex = TopoDS::Vertex(vertices->FindKey(i));
        gp_Pnt currentGeometricPoint = BRep_Tool::Pnt(vertex);
        this->pointArray->SetValue(i, currentGeometricPoint);
        this->aisShapeVertices = new AIS_Shape(vertex);
        this->aisShapeVertices->SetDisplayMode(AIS_DisplayMode::AIS_Shaded);
        this->myQtViewer->Context()->Display(this->aisShapeVertices, true);
    }
}


void MyMainWindow::onlyEdgeView()
{
    //Shaded
    this->aisShape = new AIS_Shape(*this->topoShape);
    this->aisShape->SetColor(this->defaultColor);
    this->aisShape->SetDisplayMode(AIS_DisplayMode::AIS_Shaded);

    //WireFrame
    this->aisShape = new AIS_Shape(*this->topoShape);
    this->aisShape->SetColor(Quantity_Color(0.0, 1.0, 1.0, Quantity_TOC_RGB));
    this->aisShape->SetDisplayMode(AIS_DisplayMode::AIS_WireFrame);
    this->myQtViewer->Context()->Display(this->aisShape, true);
}


void MyMainWindow::loadFile()
{
    FileDialog fileDialog;
    fileDialog.exec();
    QStringList qStrList = fileDialog.selectedFiles();
    if (qStrList.empty())
    {
        fileDialog.close();
        return;
    }
    else if (!this->topoShape->IsNull())
    {
        // reset the image when file dialog has valid input and image already exists
        // this is an implicit reset, not clicked manually by user
        this->reset(false);
    }

    std::vector<std::string> inFilePaths;
    auto numCADFiles = Utility::EnumerateCADFiles(qStrList.at(0).toStdString().c_str(), inFilePaths);
    this->myQtViewer->setFilePaths(inFilePaths);
    // visualization is applicable for single file inside a folder path or file path
    if (1 == inFilePaths.size()) // -1 == numCADFiles // QFileInfo(qStrList.at(0)).isFile()
    {
        auto singleFilePath = inFilePaths.at(0).c_str();
        std::unique_ptr<OCC_Computation> occObj = std::make_unique<OCC_Computation>(singleFilePath);
        bool status = occObj->ReadFile(singleFilePath) ? false : true;
        if (status && !occObj->GetShape()->IsNull())
        {
            this->topoShape = occObj->GetShape();
            this->aisShape = new AIS_Shape(*this->topoShape);
            this->aisShape->Color(this->defaultColor);
            this->myQtViewer->Context()->Display(this->aisShape, AIS_Shaded, 0, true);
            // scale the model as per OpenGL window size
            this->myQtViewer->View()->FitAll();
            this->myQtViewer->View()->Redraw();
            qDebug() << "\nCAD file successfully loaded!";
        }
        else
        {
            // error in file reading
            qDebug() << "\nFailed to read CAD file!";
        }
    }
    else
    {
        QMessageBox::information(this, "Selected Directory", "Visualization Disabled!\nYou Selected Directory: " + qStrList.at(0) + \
            " (" + QString::number(numCADFiles) + " Files)");
    }

    this->myQtViewer->View()->Invalidate();
    this->myQtViewer->update();
    fileDialog.close();
}


void MyMainWindow::onShowCostEstimationForm() {

    costEstimationForm = new QWidget();
    costEstimationForm->setWindowTitle("Cost Estimation");
    QSize sz(myQtViewer->window()->size().width() / 4, (myQtViewer->window()->size().height() * 0.4));
    costEstimationForm->setBaseSize(sz);
    costEstimationForm->setMaximumSize(sz);
    costEstimationForm->setMinimumSize(sz);

    // Main Layout
    QVBoxLayout* mainLayout = new QVBoxLayout(costEstimationForm);
    QGridLayout* gridLayout = new QGridLayout();  // Grid Layout for structured UI

    // Raw-Material ComboBox
    QLabel* rawMaterialLabel = new QLabel("Raw Material:");
    rawMaterialComboBox = new QComboBox();
    QStringList materialNames = { "Tough Steel", "Mild Steel", "Cast Iron", "Alloy Steel", "Stainless Steel", "Bronze", "Aluminium" };
    rawMaterialComboBox->addItems(materialNames);
    rawMaterialComboBox->setCurrentIndex(readerMachiningParameters.getMaterial()); // Default value
    gridLayout->addWidget(rawMaterialLabel, 0, 0);
    gridLayout->addWidget(rawMaterialComboBox, 0, 1);

    // Raw-Material Rate LineEdit
    QLabel* rawMaterialRateLabel = new QLabel("Raw Material Rate ₹:");
    rawMaterialRateLineEdit = new QLineEdit();
    rawMaterialRateLineEdit->setPlaceholderText(QString::number(readerMachiningParameters.getRawMaterialRate(), 'f', 2));
    rawMaterialRateLineEdit->setText(QString::number(readerMachiningParameters.getRawMaterialRate(), 'f', 2));
    gridLayout->addWidget(rawMaterialRateLabel, 1, 0);
    gridLayout->addWidget(rawMaterialRateLineEdit, 1, 1);

    // Add to main layout
    mainLayout->addLayout(gridLayout);

    // --- Map material to default rate ---
    QMap<QString, double> materialRates = {
        { "Tough Steel", 85.0 },
        { "Mild Steel", 70.0 },
        { "Cast Iron", 35.0 },
        { "Alloy Steel", 95.0 },
        { "Stainless Steel", 210.0 },
        { "Bronze", 450.0 },
        { "Aluminium", 290.0 }
    };

    // Get material name and rate from file
    int defaultMaterialIndex = readerMachiningParameters.getMaterial();
    QString defaultMaterialName = materialNames[defaultMaterialIndex];
    double defaultMaterialRate = readerMachiningParameters.getRawMaterialRate();

    // Override the map with file-loaded rate
    materialRates[defaultMaterialName] = defaultMaterialRate;

    // --- Connect ComboBox selection change to update rate ---
    connect(rawMaterialComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [=](int index) {
        QString selectedMaterial = rawMaterialComboBox->currentText();
        if (materialRates.contains(selectedMaterial)) {
            double rate = materialRates[selectedMaterial];
            rawMaterialRateLineEdit->setText(QString::number(rate, 'f', 2));
        }
        });


    // Cutter-Diameter ComboBox
    QLabel* cutterDiameterLabel = new QLabel("Cutter Diameter (mm):");
    cutterDiameterComboBox = new QComboBox();
    cutterDiameterComboBox->addItems({ "6", "10", "20", "25", "32", "40", "50", "63", "80", "100", "120", "160", "220", "250", "300" });
    cutterDiameterComboBox->setCurrentIndex(readerMachiningParameters.getCutterDiameterLegendKey());
    gridLayout->addWidget(cutterDiameterLabel, 2, 0);
    gridLayout->addWidget(cutterDiameterComboBox, 2, 1);

    // Rough-Machining-Depth Input
    QLabel* rawMachiningDepthLabel = new QLabel("Rough Machining Depth (mm):");
    roughMachiningDepthLineEdit = new QLineEdit();
    roughMachiningDepthLineEdit->setText(QString::number(readerMachiningParameters.getRoughMachiningDepth(), 'f', 2));
    gridLayout->addWidget(rawMachiningDepthLabel, 3, 0);
    gridLayout->addWidget(roughMachiningDepthLineEdit, 3, 1);


    // Rough-Machining-Depth Input
    QLabel* finishingMachiningDepthLabel = new QLabel("Finishing Machining Depth (mm):");
    finishingMachiningDepthLineEdit = new QLineEdit();
    finishingMachiningDepthLineEdit->setText(QString::number(readerMachiningParameters.getFinishingMachiningDepth(), 'f', 2));
    gridLayout->addWidget(finishingMachiningDepthLabel, 4, 0);
    gridLayout->addWidget(finishingMachiningDepthLineEdit, 4, 1);

    // Milling Rate Input
    QLabel* millingRateLabel = new QLabel("Milling Rate ₹:");
    millingRateLineEdit = new QLineEdit();
    millingRateLineEdit->setText(QString::number(readerMachiningParameters.getMillingRate(), 'f', 2));
    gridLayout->addWidget(millingRateLabel, 5, 0);
    gridLayout->addWidget(millingRateLineEdit, 5, 1);

    // Turning Rate Input
    QLabel* turningRateLabel = new QLabel("Turning Rate ₹:");
    turningRateLineEdit = new QLineEdit();
    turningRateLineEdit->setText(QString::number(readerMachiningParameters.getTurningRate(), 'f', 2));
    gridLayout->addWidget(turningRateLabel, 6, 0);
    gridLayout->addWidget(turningRateLineEdit, 6, 1);

    // Drilling Rate Input
    QLabel* drillingRateLabel = new QLabel("Drilling Rate ₹:");
    drillingRateLineEdit = new QLineEdit();
    drillingRateLineEdit->setText(QString::number(readerMachiningParameters.getDrillingRate(), 'f', 2));
    gridLayout->addWidget(drillingRateLabel, 7, 0);
    gridLayout->addWidget(drillingRateLineEdit, 7, 1);

    // Tapping Rate Input
    QLabel* tappingRateLabel = new QLabel("Tapping Rate ₹:");
    tappingRateLineEdit = new QLineEdit();
    tappingRateLineEdit->setText(QString::number(readerMachiningParameters.getTappingRate(), 'f', 2));
    gridLayout->addWidget(tappingRateLabel, 8, 0);
    gridLayout->addWidget(tappingRateLineEdit, 8, 1);

    // Buttons Layout
    QHBoxLayout* buttonLayout = new QHBoxLayout();

    // Estimate Button
    QPushButton* estimateButton = new QPushButton("Estimate Cutting Cost");
    buttonLayout->addWidget(estimateButton);

    // Save New Parameters Button
    QPushButton* saveButton = new QPushButton("Save New Parameters");
    buttonLayout->addWidget(saveButton);

    resultLabel = new QLabel(this);
    resultLabel->setText("Estimated Cost: ₹ " + QString::number(00, 'f', 2));

    // Adding Layouts to Main Layout
    mainLayout->addLayout(gridLayout);
    mainLayout->addWidget(resultLabel);
    mainLayout->addLayout(buttonLayout);
    connect(estimateButton, &QPushButton::clicked, this, [this]() { this->onEstimateButtonClicked(); });
    // Connect
    connect(estimateButton, &QPushButton::clicked, this, [this]() { this->onEstimateButtonClicked(); });
    connect(saveButton, &QPushButton::clicked, this, [this]() { this->onSaveButtonClicked(); });


    costEstimationForm->show();
}


void MyMainWindow::onEstimateButtonClicked() {
    int manufacturingCost = 0.0f;
    try
    {
        bool status = true;
        //const char* path = "D:\\Data\\Z_Feature_Recognition\\STEP";
        auto filePaths = this->myQtViewer->getFilePaths();
        if (filePaths->empty())
        {
            // error in file reading
            qDebug() << "\nFailed to load CAD file(s)!";
            return;
        }

        // save the original buffer
        std::streambuf* initBuffer = std::cout.rdbuf();
        auto numCADFiles = filePaths->size();
        std::ofstream outStream;
        if (!Utility::CONSOLE)
        {
            auto defaultOutputPath = Utility::FolderPath((*filePaths)[0].c_str());
            QString outFilePath = QFileDialog::getSaveFileName(this, "Save Compute File", defaultOutputPath.c_str(), "(*.txt)");
            if (outFilePath.isEmpty()) return;
            outStream.open(outFilePath.toStdString());
            if (outStream.good())
            {
                // redirect cout to write to the file
                std::cout.rdbuf(outStream.rdbuf());
            }
        }

        OSD_Timer timer;
        timer.Start();

        std::vector<std::array<float, 3>> slotDimensions;
        std::vector<std::pair<float, float>> holeDimensions;
        std::vector<std::pair<float, float>> handleDimensions;
        std::array<float, 3> boxDimensions{ 0.0f };
        std::vector<std::pair<float, float>> threadedHoleDimensions;




        if (1 == numCADFiles)
        {
            // file path or folder path with single file
            std::unique_ptr<AS_FeatureRecognition> asObj;
            auto fileExtn = Utility::FileExtn(filePaths->at(0).c_str());
            if ("step" == fileExtn || "stp" == fileExtn)
            {
                // for single CAD file, TopoDS_Shape is already populated during loading operation
                // use this shape only for STEP files because AS_FeatureRecognition::IGESReader()
                // is different from OCC_Computation::ReadFile() for IGES files
                // TopoGeomMaps are different in both the class methods and result in downstream errors
                // if the same shape is used
                assert(!topoShape->IsNull());
                asObj = std::make_unique<AS_FeatureRecognition>(filePaths->at(0).c_str(), topoShape.get());
            }
            else if ("iges" == fileExtn || "igs" == fileExtn)
                asObj = std::make_unique<AS_FeatureRecognition>(filePaths->at(0).c_str());
#if OCC
            auto occObj = asObj->GetCADPart();
            bool status = occObj->ReadFile() ? false : true;
            if (status) status &= occObj->GenerateTopologyMaps();
            //if (status) status &= occObj->WriteFile(false) ? false : true;
            if (status) status &= occObj->Compute();
#else
            status = asObj->RunFeatureRecognition(filePaths->at(0),
                slotDimensions, holeDimensions, threadedHoleDimensions, handleDimensions, boxDimensions);
            std::unique_ptr<CostEstimation> ceObj = std::make_unique<CostEstimation>
                (slotDimensions, holeDimensions, handleDimensions, boxDimensions, readerMachiningParameters.getMaterial(), readerMachiningParameters.getRoughMachiningDepth(),
                    readerMachiningParameters.getFinishingMachiningDepth(), readerMachiningParameters.getCutterDiameter(), readerMachiningParameters.getRawMaterialRate(), readerMachiningParameters.getMillingRate(),
                    readerMachiningParameters.getTurningRate(), readerMachiningParameters.getDrillingRate(), readerMachiningParameters.getTappingRate());
            manufacturingCost = ceObj->RunCostEstimation(&slotDimensions, &holeDimensions, &handleDimensions, &boxDimensions);
#endif
        }
        else if (numCADFiles > 0)
        {
            // folder path
            // std::make_unique supports default initialization
            // new[] supports parameterized initialization
            std::unique_ptr<AS_FeatureRecognition[]> asObjs = std::make_unique<AS_FeatureRecognition[]>(numCADFiles); //readerMachiningParameters

            //rawMaterial, roughCuttingDepth, finishCuttingDepth,
            //cutterDiameter, rawMaterialRate, millingRate, turningRate, drillingRate, tappingRate
            std::unique_ptr<CostEstimation[]> ceObjs(new CostEstimation[numCADFiles]{ CostEstimation(readerMachiningParameters.getMaterial(),readerMachiningParameters.getRoughMachiningDepth(),
                readerMachiningParameters.getFinishingMachiningDepth(),readerMachiningParameters.getCutterDiameter(),readerMachiningParameters.getMillingRate(),
                readerMachiningParameters.getTurningRate(),readerMachiningParameters.getDrillingRate(), readerMachiningParameters.getTappingRate()) });
            for (int16_t i = 0; i < numCADFiles; ++i)
            {
                status &= asObjs[i].RunFeatureRecognition(filePaths->at(i),
                    slotDimensions, holeDimensions, threadedHoleDimensions, handleDimensions, boxDimensions);
                if (status)
                    manufacturingCost += ceObjs[i].RunCostEstimation
                    (&slotDimensions, &holeDimensions, &handleDimensions, &boxDimensions);
            }

            if (status)
                qDebug() << "\nCAD files successfully loaded!";
        }
        else
            status &= false;
        std::cout << "\n\n";
        timer.Stop();
        timer.Show();

        // reset the buffer to print to the console again
        if (!Utility::CONSOLE)
        {
            std::cout.rdbuf(initBuffer);
            outStream.close();
        }
        if (status)
            qDebug() << "\nCAD manufacturing cost estimated successfully!";
    }
    catch (...)
    {
        qDebug() << "\nEXCEPTION-CAUGHT !!!";
    }


    // Get the values selected by the user
    QString rawMaterial = rawMaterialComboBox->currentText();
    QString cutterDiameterStr = cutterDiameterComboBox->currentText();
    float rawMachiningDepth = roughMachiningDepthLineEdit->text().toFloat();

    if (rawMachiningDepth <= 0)
        QMessageBox::warning(this, "Invalid Input", "Please enter a valid machining depth.");


    // CostEstimation estimateMachingCost;
     //QString tcs = QString::number(manufacturingCost, 'f', 2); // format to 2 decimal places
     //resultLabel->setText("Estimated Cost: ₹ " + tcs);

    QString styledText = "<span style='font-size:12pt; color:#d32f2f; font-weight:bold;'>"
        "Estimated Cost: ₹ " + QString::number(manufacturingCost, 'f', 2) +
        "</span>";

    resultLabel->setText(styledText);

    // costEstimationForm->close();

}


void MyMainWindow::onSaveButtonClicked() {
    // Update parameters from user input
    readerMachiningParameters.setMaterial(rawMaterialComboBox->currentIndex());
    readerMachiningParameters.setRawMaterialRate(rawMaterialRateLineEdit->text().toDouble());
    readerMachiningParameters.setCutterDiameter(cutterDiameterComboBox->currentText().toDouble());
    readerMachiningParameters.setRoughMachiningDepth(roughMachiningDepthLineEdit->text().toDouble());
    readerMachiningParameters.setFinishingMachiningDepth(finishingMachiningDepthLineEdit->text().toDouble());
    readerMachiningParameters.setMillingRate(millingRateLineEdit->text().toDouble());
    readerMachiningParameters.setTurningRate(turningRateLineEdit->text().toDouble());
    readerMachiningParameters.setDrillingRate(drillingRateLineEdit->text().toDouble());
    readerMachiningParameters.setTappingRate(tappingRateLineEdit->text().toDouble());

    // Save to file
    readerMachiningParameters.saveMachingParameters();

    // Show confirmation message
    QMessageBox::information(costEstimationForm, "Success", "Machining parameters saved successfully!");

}


void MyMainWindow::compute(const std::shared_ptr<TopoDS_Shape>& topoShape)
{
    auto filePaths = this->myQtViewer->getFilePaths();
    if (filePaths->empty())
    {
        // error in file reading
        qDebug() << "\nFailed to load CAD file(s)!";
        return;
    }

    else
        //costestimatation form contains machining parameters
        onShowCostEstimationForm();
}


void MyMainWindow::reset(bool flag)
{
    if (!this->aisShape.IsNull() && !this->topoShape->IsNull())
    {
        this->wireframeCheckBox.setCheckState(Qt::CheckState::Unchecked);
        this->verticesCheckBox.setCheckState(Qt::CheckState::Unchecked);
        this->edgeCheckBox.setCheckState(Qt::CheckState::Unchecked);
        this->myQtViewer->Context()->Remove(this->aisShape, true);
        this->myQtViewer->Context()->Remove(this->aisShapeVertices, true);
        this->topoShape->Nullify();
        this->aisShape.Nullify();
        this->aisShapeVertices.Nullify();
        this->vertices->Clear();
        this->myQtViewer->setFilePaths({});
        this->myQtViewer->View()->Invalidate();
        this->myQtViewer->update();
        qDebug() << "\nResetting CAD file(s)!";
    }
    else if (!this->myQtViewer->getFilePaths()->empty())
    {
        this->myQtViewer->setFilePaths({});
        qDebug() << "\nResetting CAD file(s)!";
    }

    // reset the radio button's state only for explicit reset button click
    if (flag && 0 != this->myButtonGroup.checkedButton())
    {
        this->myButtonGroup.button(2)->setChecked(true);
        Utility::CONSOLE = false;
    }
}


void MyMainWindow::checkAndView()
{
    auto isWireframeModeChecked = (this->wireframeCheckBox.checkState() == Qt::CheckState::Checked);
    auto isVerticesModeChecked = (this->verticesCheckBox.checkState() == Qt::CheckState::Checked);
    auto isEdgeModeChecked = (this->edgeCheckBox.checkState() == Qt::CheckState::Checked);

    try
    {
        if (!this->aisShape.IsNull() && !this->topoShape->IsNull())
        {
            if (isWireframeModeChecked && isEdgeModeChecked && isVerticesModeChecked)
            {
                this->firstLoadView();
                this->myQtViewer->Context()->Display(this->aisShape, AIS_WireFrame, 0, true);
                this->onlyVerticesView();
                this->onlyEdgeView();
            }
            else if (isWireframeModeChecked && isEdgeModeChecked && !isVerticesModeChecked)
            {
                this->firstLoadView();
                this->myQtViewer->Context()->Display(this->aisShape, AIS_WireFrame, 0, true);
                this->onlyEdgeView();
            }
            else if (isWireframeModeChecked && !isEdgeModeChecked && isVerticesModeChecked)
            {
                this->firstLoadView();
                this->myQtViewer->Context()->Display(this->aisShape, AIS_WireFrame, 0, true);
                this->onlyVerticesView();
            }
            else if (!isWireframeModeChecked && isEdgeModeChecked && isVerticesModeChecked)
            {
                this->firstLoadView();
                this->onlyVerticesView();
                this->onlyEdgeView();

            }
            else if (isWireframeModeChecked && !isEdgeModeChecked && !isVerticesModeChecked)
            {
                this->firstLoadView();
                this->myQtViewer->Context()->Display(this->aisShape, AIS_WireFrame, 0, true);
            }
            else if (!isWireframeModeChecked && isEdgeModeChecked && !isVerticesModeChecked)
            {
                this->firstLoadView();
                this->onlyEdgeView();
            }
            else if (!isWireframeModeChecked && !isEdgeModeChecked && isVerticesModeChecked)
            {
                this->firstLoadView();
                this->onlyVerticesView();
            }
            else
            {
                this->firstLoadView();
            }
        }
        else
        {
            if (isWireframeModeChecked || isEdgeModeChecked || isVerticesModeChecked)
            {
                this->wireframeCheckBox.setCheckState(Qt::CheckState::Unchecked);
                this->verticesCheckBox.setCheckState(Qt::CheckState::Unchecked);
                this->edgeCheckBox.setCheckState(Qt::CheckState::Unchecked);
                qDebug() << "\nFailed to load CAD file!";
            }
        }
    }
    catch (...)
    {
        qDebug() << "\nEXCEPTION-CAUGHT !!!";
    }
}


void MyMainWindow::closeEvent(QCloseEvent* event)
{

    QSettings settings("YourCompany", "YourApp");
    settings.setValue("printOption", this->myButtonGroup.checkedId());

    QMainWindow::closeEvent(event);

    // Close the Cost Estimation form when the parent window is closed
    if (costEstimationForm) costEstimationForm->close();  // Close the child form
    else if (!costEstimationForm)event->accept();
    else event->accept();
}