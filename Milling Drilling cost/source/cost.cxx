#include <algorithm>
// project
#include <cost.hxx>

constexpr int16_t MillingFeed[2][7][15] =
{
    // rough - 15 columns
    {
        { 25, 20, 31, 24, 76, 76, 61, 73, 57, 61, 64, 67, 63, 61, 61 },
        { 53, 42, 63, 50, 153, 153, 122, 146, 115, 122, 127, 134, 125, 122, 122 },
        { 33, 27, 40, 32, 172, 172, 138, 164, 129, 138, 143, 150, 141, 138, 138 },
        { 45, 36, 54, 43, 153, 153, 122, 146, 115, 122, 127, 134, 125, 122, 122 },
        { 49, 39, 59, 47, 134, 134, 107, 127, 100, 107, 111, 117, 109, 107, 107 },
        { 56, 45, 67, 54, 134, 134, 107, 127, 100, 107, 111, 117, 109, 107, 107 },
        { 96, 76, 115, 92, 210, 210, 168, 200, 158, 168, 175, 184, 172, 168, 168 }
    },
    // finish - 15 columns
    {
        { 19, 15, 23, 18, 57, 57, 46, 55, 43, 46, 48, 50, 47, 46, 46 },
        { 45, 36, 54, 43, 115, 115, 92, 109, 86, 92, 96, 100, 94, 92, 92 },
        { 29, 23, 34, 28, 134, 134, 107, 127, 100, 107, 111, 117, 109, 107, 107 },
        { 35, 28, 42, 34, 115, 115, 92, 109, 86, 92, 96, 100, 94, 92, 92 },
        { 35, 28, 42, 34, 105, 105, 84, 100, 79, 84, 88, 92, 86, 84, 84 },
        { 40, 32, 48, 38, 105, 105, 84, 100, 79, 84, 88, 92, 86, 84, 84 },
        { 64, 51, 76, 61, 162, 162, 130, 155, 122, 130, 135, 142, 133, 130, 130 }
    }
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

constexpr int16_t TurningFeed[2][7][12] =
{
    // rough - 12 columns
    {
        {25, 13,  8,  6,  5,  4,  4,  3,  3,  3,  2,  2},
        {51, 25, 17, 13, 10,  8,  7,  6,  6,  5,  5,  4},
        {57, 29, 19, 14, 11, 10,  8,  7,  6,  6,  5,  5},
        {51, 25, 17, 13, 10,  8,  7,  6,  6,  5,  5,  4},
        {45, 22, 15, 11,  9,  7,  6,  6,  5,  4,  4,  4},
        {45, 22, 15, 11,  9,  7,  6,  6,  5,  4,  4,  4},
        {70, 35, 23, 18, 14, 12, 10,  9,  8,  7,  6,  6}
    },
    // finish - 12 columns
    {
        {19, 10,  6,  5,  4,  3,  3,  2,  2,  2,  2,  2},
        {38, 19, 13, 10,  8,  6,  5,  5,  4,  4,  3,  3},
        {45, 22, 15, 11,  9,  7,  6,  6,  5,  4,  4,  4},
        {38, 19, 13, 10,  8,  6,  5,  5,  4,  4,  3,  3},
        {35, 18, 12,  9,  7,  6,  5,  4,  4,  4,  3,  3},
        {35, 18, 12,  9,  7,  6,  5,  4,  4,  4,  3,  3},
        {54, 27, 18, 14, 11,  9,  8,  7,  6,  5,  5,  5}
    }
};

const std::map<float, float> DiaToDrillDepth = {
    {1,24},{1.5,30},{1.6,32.5},{1.8,34.22222222},{2,36},
    {2.1,35.42857143},{2.2,39.27272727},{2.3,38.73913043},{2.4,42.5},{2.5,42},
    {2.6,41.53846154},{2.7,45.22222222},{2.78,44.8705036},{2.8,44.78571429},{2.9,44.37931034},{3,44},
    {3.1,47.61290323},{3.2,47.25},{3.3,46.90909091},{3.4,50.47058824},{3.5,50.14285714},
    {3.6,49.83333333},{3.7,49.54054054},{3.8,54.31578947},{3.9,54.02564103},{4,53.75},
    {4.1,53.48780488},{4.2,53.23809524},{4.3,57.93023256},{4.37,57.75514874},{4.4,57.68181818},{4.5,57.44444444},
    {4.6,57.2173913},{4.7,57},{4.75,62.94736842},{4.8,62.83333333},{4.9,62.6122449},{5,62.4},
    {5.8,66.82758621},{5.9,66.66101695},{6,66.5},
    {6.1,73.32786885},{6.2,73.16129032},{6.3,73},{6.4,72.84375},{6.5,72.69230769},
    {6.6,72.54545455},{6.7,72.40298507},{6.75,79.22222222},{6.8,79.14705882},{6.9,79},{7,78.85714286},
    {7.1,78.71830986},{7.2,78.58333333},{7.3,78.45205479},{7.4,78.32432432},{7.5,78.2},
    {7.6,84.86842105},{7.7,84.74025974},{7.8,84.61538462},{7.9,84.49367089},{8,84.375},
    {8.1,84.25925926},{8.2,84.14634146},{8.3,84.03614458},{8.4,83.92857143},{8.5,83.82352941},
    {8.6,90.41860465},{8.7,90.31034483},{8.75,90.25714286},{8.8,90.20454545},{8.9,90.1011236},{9,90},
    {9.1,89.9010989},{9.2,89.80434783},{9.3,89.70967742},
    {9.4,89.61702128},{9.5,89.52631579},{9.6,96.0625},{9.7,95.96907216},{9.8,95.87755102},{9.9,95.78787879},{10,95.7},
    {10.2,95.52941176},{10.25,95.48780488},{10.5,95.28571429},{10.8,102.7037037},{11,102.5454545},
    {11.1,102.4684685},{11.5,102.173913},{12,109.4166667},
    {12.5,109.08},{12.7,108.9527559},{13,108.7692308}
};

bool isDrillDepthLessHole(float holeRadius, float holeHeight) {
    for (const auto& [dia, depth] : DiaToDrillDepth) {
        if (holeRadius <= dia && std::floor(holeRadius) == std::floor(dia)) {
            std::cout << "  Drill Radius:" << holeRadius << " Drill Height:" << holeHeight << std::endl;
            return holeHeight <= depth;
        }
    }
    return false; // If no matching diameter found
}



// TOUGH_STEEL=7.85, MILD_STEEL=7.86, CAST_IRON=7.2, ALLOY_STEEL=7.8, STAINLESS_STEEL=7.95, BRONZE=8.73, ALUMINIUM=2.7
constexpr float MaterialDensity[7] = { 7.85, 7.86, 7.2, 7.8, 7.95, 8.73, 2.7 };


CostEstimation::CostEstimation(
    //std::shared_ptr<OcctQtViewer> myQtViewer,
    int16_t rawMaterial,
    float roughCuttingDepth, float finishCuttingDepth,
    int16_t cutterDiameter, int16_t rawMaterialRate,
    int16_t millingRate, int16_t turningRate,
    int16_t drillingRate, int16_t tappingRate)
    : CostEstimation(
       // myQtViewer,
        {}, {}, {}, {},
        rawMaterial, roughCuttingDepth, finishCuttingDepth,
        cutterDiameter, rawMaterialRate, millingRate, turningRate, drillingRate, tappingRate) {}


CostEstimation::CostEstimation(
   // std::shared_ptr<OcctQtViewer> myQtViewer,
    const std::vector<std::array<float, 3>>& slotDimensions,
    const std::vector<std::pair<float, float>>& holeDimensions,
    const std::vector<std::pair<float, float>>& handleDimensions,
    const std::array<float, 3>& boxDimensions,
    int16_t rawMaterial,
    float roughCuttingDepth, float finishCuttingDepth,
    int16_t cutterDiameter, int16_t rawMaterialRate,
    int16_t millingRate, int16_t turningRate,
    int16_t drillingRate, int16_t tappingRate)
{
   // this->myQtViewer = myQtViewer;
    this->m_Handles = std::make_shared <std::vector<std::pair<float, float>>>(handleDimensions);
    this->m_Holes = std::make_shared <std::vector<std::pair<float, float>>>(holeDimensions);
    this->m_Slots = std::make_shared<std::vector<std::array<float, 3>>>(slotDimensions);
    if (!boxDimensions.empty())
        this->m_BoxVolume = boxDimensions[0] * boxDimensions[1] * boxDimensions[2];
    this->m_RawMaterial = rawMaterial;
    this->m_CutterDiameter = cutterDiameter;
    this->m_RoughMillingCuttingDepth = roughCuttingDepth;
    this->m_FinishMillingCuttingDepth = finishCuttingDepth;
    // map rawMaterial string to RawMaterial enum
    // map cutterDiameter short int to CutterDiameter enum
    int material = rawMaterial; //std::underlying_type<CostEstimation::RawMaterial>::type(RawMaterial::ALLOY_STEEL);
    int16_t diameter = cutterDiameter;// std::underlying_type<CostEstimation::CutterDiameter>::type(CutterDiameter::EIGHTY);

    // Iterate through the map to find the key for the cutter Diameter value
    for (const auto& pair : CutterDiameterLegend) {
        if (pair.second == cutterDiameter) {
            diameter = pair.first;
            break;
        }
    }


    this->m_RoughCuttingFeed = MillingFeed[0][material][diameter];
    this->m_FinishCuttingFeed = MillingFeed[1][material][diameter];
    this->m_MaterialDensity = MaterialDensity[material];
    this->m_MaterialRate = rawMaterialRate;
    this->m_MillingRate = millingRate;
    this->m_TurningRate = turningRate;
    this->m_DrillingRate = drillingRate;
    this->m_TappingRate = tappingRate;
}


CostEstimation::~CostEstimation() {}

float CostEstimation::RawMaterialCost(std::array<float, 3>* boxDimensions)
{
    float totalRawMaterialRate = 0.0f;
    auto weight = boxDimensions->at(0) * boxDimensions->at(1) * boxDimensions->at(2) * this->m_MaterialDensity;
    weight = weight / 1000000;
    totalRawMaterialRate = weight * this->m_MaterialRate;

    std::cout << std::endl;
    std::cout << "CostEstimation: \n" << std::endl;
    std::cout << " Material: " << materialLegend.at(this->m_RawMaterial) << std::endl;
    std::cout << " Total Raw Material Rate: " << totalRawMaterialRate << " Rs" << std::endl;

    return totalRawMaterialRate;
}

float CostEstimation::MillingCost(std::vector<std::array<float, 3>>* slotDimensions)
{
    float totalMillingCost = 0.0f; int slotcnt = 0;



    if (slotDimensions && slotDimensions->empty())return 0;

    // --- Milling Cost --- 
    std::cout << "\n Milling Cost: " << std::endl;
    for (const auto& slot : *slotDimensions)
    {
        std::cout << "  Slot:" << slotcnt << std::endl;
        auto minLength = std::min(slot.at(0), std::min(slot.at(1), slot.at(2)));
        auto maxLength = std::max(slot.at(0), std::max(slot.at(1), slot.at(2)));

        // length, width, height of slot
        float length{ 0 }, width{ 0 }, height{ 0 };
        length = maxLength + this->m_CutterDiameter;
        for (auto side : slot) {
            if (!(std::fabs(side - maxLength) < 1e-6) && !(std::fabs(side - minLength) < 1e-6))width = side;
        }
        height = minLength;
        float depthOfMilling = minLength;
        std::cout << "  Length: " << length << ", Width: " << width << ", Depth Of Milling: " << depthOfMilling << std::endl;

        auto depthofcutforFinishMilling = this->m_FinishMillingCuttingDepth;
        auto countoftravelforfinishMilling = 2;
        auto totalDepthApplicableForFinishMilling = countoftravelforfinishMilling * depthofcutforFinishMilling;
        auto totalDepthApplicableForRawMilling = depthOfMilling - totalDepthApplicableForFinishMilling;
        std::cout << "  Depth of cut for Finish Milling: " << depthofcutforFinishMilling << ", Count of travel for Finish Milling: " << countoftravelforfinishMilling << std::endl;
        std::cout << "  Total Depth Applicable For Finish Milling: " << totalDepthApplicableForFinishMilling << ", Total Depth Applicable For Raw Milling: " << totalDepthApplicableForRawMilling << std::endl;

        // rough milling
        auto noOfCutsRoughMilling = (width / this->m_CutterDiameter) * (totalDepthApplicableForRawMilling / this->m_RoughMillingCuttingDepth);
        auto timeForRoughMilling = (length / this->m_RoughCuttingFeed) * noOfCutsRoughMilling;
        std::cout << "  No Of Cuts Rough Milling: " << noOfCutsRoughMilling << ", Rough Cutting Feed:" << this->m_RoughCuttingFeed << "\n  Time For Rough Milling: " << timeForRoughMilling << std::endl;

        // finish milling
        auto noOfCutsFinishMilling = (width / this->m_CutterDiameter) * (totalDepthApplicableForFinishMilling / this->m_FinishMillingCuttingDepth);
        auto timeForFinishMilling = (length / this->m_FinishCuttingFeed) * noOfCutsFinishMilling;
        std::cout << "  No Of Cuts Finish Milling: " << noOfCutsFinishMilling << ", Finish Cutting Feed:" << this->m_FinishCuttingFeed << "\n  Time For Finish Milling: " << timeForFinishMilling << std::endl;

        // total time for milling
        auto totalTimeForMillingMin = timeForRoughMilling + timeForFinishMilling;
        auto totalTimeForMillingHrs = totalTimeForMillingMin / 60;
        std::cout << "\n  Total Time For Milling (in min.):" << totalTimeForMillingMin << ", Milling Rate (per hr):" << this->m_MillingRate << " Rs" << std::endl;

        // total cost for milling slot
        auto totalCostForSlotMilling = totalTimeForMillingHrs * this->m_MillingRate;

        totalMillingCost = totalMillingCost + totalCostForSlotMilling;
        std::cout << "\n  Total Milling Cost for Slot " << slotcnt << ": " << totalMillingCost << " Rs\n" << std::endl;

        ++slotcnt;
    }

    return totalMillingCost;
}

float CostEstimation::DrillingTappingCost(std::vector<std::pair<float, float>>* holeDimensions)
{
    static bool alreadyAsked = false;
    float totalDrillingCost = 0.0f;
    float totalTappingCost = 0.0f;

    std::cout << "\nDrilling Cost:" << std::endl;
    std::cout << "\nDrills: " << holeDimensions->size() << std::endl;

    noDrillTooShallow = 0;
    if (!holeDimensions->empty()) {
        for (const auto& [radius, height] : *holeDimensions) {
            if (!isDrillDepthLessHole(radius, height)) {
                ++noDrillTooShallow;
            }
        }
    }

    int holeCount = static_cast<int>(holeDimensions->size());
    float extraCostDrill = this->m_DrillingRate / 2.0f;


    

    std::cout << "[DEBUG] DrillingTappingCost() called" << std::endl;

    if (alreadyAsked) {
        std::cout << "[DEBUG] Already asked about Special Tool. Skipping UI." << std::endl;
    }
     useSpecialTool = false;
     float spacialToolCost = 200.0f; // Example special tool cost
    if (!alreadyAsked && noDrillTooShallow > 0) {
        alreadyAsked = true;
        // Create Form
        drillCostForm = new QWidget();
        drillCostForm->setWindowTitle("Drill Cost");
        QSize sz(300, 120);
        drillCostForm->setBaseSize(sz);
        drillCostForm->setMaximumSize(sz);
        drillCostForm->setMinimumSize(sz);

        // Layouts
        QVBoxLayout* mainLayout = new QVBoxLayout(drillCostForm);
        QLabel* decisionLabel = new QLabel("Use Special Tool for Drilling?");
        mainLayout->addWidget(decisionLabel, 0, Qt::AlignCenter);

                // Show extra cost and special tool cost
        QLabel* extraCostLabel = new QLabel(QString("Extra cost for drill:₹ %1 ").arg(extraCostDrill));
        mainLayout->addWidget(extraCostLabel, 0, Qt::AlignCenter);

        
        QLabel* spacialToolLabel = new QLabel(QString("Special Tool Cost:₹ %1 ").arg(spacialToolCost));
        mainLayout->addWidget(spacialToolLabel, 0, Qt::AlignCenter);


        QHBoxLayout* buttonLayout = new QHBoxLayout();

        QPushButton* spacialToolYButton = new QPushButton("Yes - Use Special Tool");
        QPushButton* spacialToolNButton = new QPushButton("No - Use Default Flow");

        buttonLayout->addWidget(spacialToolYButton);
        buttonLayout->addWidget(spacialToolNButton);
        mainLayout->addLayout(buttonLayout);
        
        // Event Loop to wait for user's choice
        QEventLoop loop;

        connect(spacialToolYButton, &QPushButton::clicked, this, [&]() {
            std::cout << "User selected YES" << std::endl;
            useSpecialTool = true;
            //alreadyAsked = false;
            drillCostForm->close();
            
            loop.quit();
            });

        connect(spacialToolNButton, &QPushButton::clicked, this, [&]() {
            std::cout << "User selected NO" << std::endl;
            useSpecialTool = false;
           // alreadyAsked = false;
            drillCostForm->close();
            loop.quit();
            });

        drillCostForm->show();
        loop.exec();  // Pause until button is clicked
    }

    // Apply modified cost if special tool is selected
    if (useSpecialTool) {
        extraCostDrill *= spacialToolCost; // Example: increase cost more if using special tool
    }

    totalDrillingCost = (this->m_DrillingRate * holeCount) + (extraCostDrill * noDrillTooShallow);

    std::cout << "Total Drilling Cost: " << totalDrillingCost << " Rs"
        << " | Drilling Rate: " << this->m_DrillingRate << " Rs"
        << " | Shallow Holes: " << noDrillTooShallow
        << " | Extra Cost per Shallow Hole: " << extraCostDrill
        << std::endl;

    return totalDrillingCost;
}






float CostEstimation::TurningCost()
{
    return 0.0f;
}


int CostEstimation::RunCostEstimation(std::vector<std::array<float, 3>>* slotDimensions,
    std::vector<std::pair<float, float>>* holeDimensions,
    std::vector<std::pair<float, float>>* handleDimensions,
    std::array<float, 3>* boxDimensions)
{
    int manufacturingCost = 0.0f;  float totalTappingCost = 0.0f;

    // --- Raw Material Cost ---
    float totalRawMaterialRate = RawMaterialCost(boxDimensions);


    // --- Milling Cost --- 
    float totalMillingCost = MillingCost(slotDimensions);


    // --- Drilling Cost ---
    float totalDrillingCost = DrillingTappingCost(holeDimensions);

    // --- Turning Cost ---
    float totalTurningCost = TurningCost();

    // ---  Total Manufacturing Cost ---
    manufacturingCost = totalRawMaterialRate + totalMillingCost + totalDrillingCost;
    std::cout << " Total Manufacturing Cost:" << manufacturingCost << " Rs" << std::endl;
    std::cout << std::endl;

    return manufacturingCost;
}


void CostEstimation::closeEvent(QCloseEvent* event)
{
    // Close the Cost Estimation form when the parent window is closed
    if (drillCostForm) drillCostForm->close();  // Close the child form
    else if (!drillCostForm)event->accept();
    else event->accept();
}



