#include <asi.hxx>
#include <common.hxx>
#include <occ.hxx>

#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <TopExp_Explorer.hxx>
#include <TopOpeBRepBuild_Tools.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>

#include <asiAlgo.h>
#include <asiAlgo_AAGIterator.h>
#include <asiAlgo_CheckDihedralAngle.h>
#include <asiAlgo_IGES.h>
#include <asiAlgo_CheckValidity.h>
#include <asiAlgo_ConvertCanonical.h>
#include <asiAlgo_ConvertCanonicalSummary.h>
#include <asiAlgo_ReadSTEPWithMeta.h>
#include <asiAlgo_RecognizeDrillHoles.h>
#include <asiAlgo_RecognizeShafts.h>
#include <asiAlgo_STEP.h>
#include <asiAlgo_Utils.h>
#include <asiData.h>
#include <asiEngine.h>
#include <ActData.h>

#include <array>
#include <set>


AS_FeatureRecognition::AS_FeatureRecognition(const char* cadFilePath)
{
    this->m_CADPart = std::make_shared<OCC_Computation>(cadFilePath);
    this->m_ConvexEdges = std::make_shared<std::map<int, TopoDS_Edge>>();
    this->m_ConcaveEdges = std::make_shared<std::map<int, TopoDS_Edge>>();
    this->m_SmoothC1Edges = std::make_shared<std::map<int, TopoDS_Edge>>();
    this->m_HoleFaces = std::make_shared<TColStd_PackedMapOfInteger>();
    //std::cout << "Constructor!" << std::endl;
}


AS_FeatureRecognition::~AS_FeatureRecognition()
{
    //std::cout << "Destructor!" << std::endl;
}


bool AS_FeatureRecognition::STEPReader()
{
    Handle(asiAlgo_STEP) stp = new asiAlgo_STEP;
    return stp->Read(TCollection_AsciiString(this->m_CADPart->GetFilePath().c_str()), false, *this->m_CADPart->GetShape().get());
    //Handle(asiAlgo_ReadSTEPWithMeta) reader = new asiAlgo_ReadSTEPWithMeta;
    //bool status = reader->Perform(this->m_FilePath.c_str());
}


bool AS_FeatureRecognition::IGESReader()
{
    Handle(asiAlgo_IGES) igs = new asiAlgo_IGES(ActAPI_ProgressEntry());
    bool success = igs->Read(TCollection_AsciiString(this->m_CADPart->GetFilePath().c_str()), *this->m_CADPart->GetShape().get());
    if (success)
    {
        success &= this->m_CADPart->GenerateTopologyGeometryMaps();
        success &= this->m_CADPart->ConvertIGESToSolid();
        //this->m_CADPart->ResetTopologyMaps();
        //success &= this->m_CADPart->GenerateTopologyGeometryMaps();
        //this->m_CADPart->ComputePartGeometry(true);
        success &= this->ConvertToCanonicalGeometry();
        // necessary because conversion of surfaces to solids changes both topology & geometry
        this->m_CADPart->ResetTopologyGeometryMaps();
        success &= this->m_CADPart->GenerateTopologyGeometryMaps();
        //this->m_CADPart->ComputePartGeometry(true);
    }
    return success;
}


bool AS_FeatureRecognition::GenerateAAG()
{
    bool success = true;
    if (this->m_CADPart->GetShape()->IsNull())
    {
        auto fileExtn = Utility::FileExtn(this->m_CADPart->GetFilePath().c_str());
        if ("step" == fileExtn || "stp" == fileExtn)
            success = this->STEPReader();
        else if ("iges" == fileExtn || "igs" == fileExtn)
            success = this->IGESReader();
    }
    if (success)
    {
        // attribute adjacency graph
        this->m_AAG = new asiAlgo_AAG(*this->m_CADPart->GetShape());
#if 0
        this->m_AAG->Dump(std::cout);
#endif
    }

    return success;
}


bool AS_FeatureRecognition::Initialize()
{
    bool status = this->GenerateAAG();
    if (!this->m_CADPart->GetTopoGeomMapsGenerated())
        status &= this->m_CADPart->GenerateTopologyGeometryMaps();
    //status &= this->m_CADPart->ComputePartGeometryMaps();
    status &= this->ClassifyEdges();
    status &= this->FindHoleFaces();
    //...
    return status;
}


bool AS_FeatureRecognition::ConvertToCanonicalGeometry()
{
    auto cadPart = this->m_CADPart;
    if (cadPart->GetShape()->IsNull()) return false;

    asiAlgo_ConvertCanonical converter;
    TopoDS_Shape modifiedShape = converter.Perform(*cadPart->GetShape(), cadPart->GetTolerance());
    asiAlgo_CheckValidity checker;
    bool status = checker.CheckBasic(modifiedShape);
    if (!status) return status;

    // feature recognition section of the code is not supposed to modify the original CAD shape
    // but in case of IGES files, feature recognition is easier with canonical (analytical) shapes
    this->m_CADPart->SetShape(std::make_shared<TopoDS_Shape>(modifiedShape));
    asiAlgo_ConvertCanonicalSummary summary = converter.GetSummary();
    summary.isValid = true;
    return status;
}


AS_FeatureRecognition::PartCategory AS_FeatureRecognition::IdentifyPartCategory()
{
    // this criteria will change as more CAD part categories are introduced
    AS_FeatureRecognition::PartCategory enumObj = AS_FeatureRecognition::PartCategory::UNDEFINED;
    for (const auto& itr : *this->m_ConcaveEdges)
    {
        BRepAdaptor_Curve curve(itr.second);
        if (GeomAbs_CurveType::GeomAbs_Line == curve.GetType())
        {
            // defines a slot
            enumObj = AS_FeatureRecognition::PartCategory::MILLING;
            break;
        }
        else if (GeomAbs_CurveType::GeomAbs_Circle == curve.GetType())
            enumObj = AS_FeatureRecognition::PartCategory::TURNING;
    }

    if (AS_FeatureRecognition::PartCategory::UNDEFINED == enumObj)
    {
        int16_t numCylindricalFaces = 0;
        std::map<int16_t, GeomAbs_CurveType> etMap;
        std::map<int16_t, GeomAbs_SurfaceType> ftMap;
        this->m_CADPart->QueryGeometryMaps(etMap, ftMap);

        for (const auto& itr : ftMap)
        {
            if (GeomAbs_SurfaceType::GeomAbs_Cylinder == itr.second)
            {
                // find only those cylindrical faces which are not categorized as holes
                if (this->m_HoleFaces->Contains(itr.first))
                    continue;
                ++numCylindricalFaces;
            }
        }
        if (numCylindricalFaces)
            enumObj = AS_FeatureRecognition::PartCategory::TURNING;
    }

    return enumObj;
}


double AS_FeatureRecognition::ComputeDihedralAngle( const TopoDS_Face& face1,
                                                    const TopoDS_Face& face2,
                                                    bool allowC1Edges,
                                                    double angularTolerance,
                                                    std::map<int, TopoDS_Edge>* convexEdges,
                                                    std::map<int, TopoDS_Edge>* concaveEdges,
                                                    std::map<int, TopoDS_Edge>* smoothC1Edges,
                                                    std::map<int, TopoDS_Edge>* undefinedEdges  )
{
    auto InsertEdges = [](  const TopTools_IndexedMapOfShape& sharedEdges,
                            const std::vector<int16_t>& edgeIDs,
                            std::map<int, TopoDS_Edge>* edgeType) -> void
    {
        for (int16_t idx = 1; idx <= sharedEdges.Extent(); ++idx)
            edgeType->insert(std::make_pair(edgeIDs.at(idx - 1), TopoDS::Edge(sharedEdges.FindKey(idx))));
    };

    TopTools_IndexedMapOfShape vMap, eMap, fMap;
    TopTools_IndexedDataMapOfShapeListOfShape veMap, efMap;
    bool status = this->m_CADPart->QueryTopologyMaps(vMap, eMap, fMap, veMap, efMap);
    if (!status) return std::numeric_limits<double>::max();

    TopTools_IndexedMapOfShape sharedEdges;
    asiAlgo_CheckDihedralAngle dihedralAngle;
    double angle = 0.0;
    asiAlgo_FeatureAngleType angleType = 
        dihedralAngle.AngleBetweenFaces(face1, face2, allowC1Edges, angularTolerance,
                                        sharedEdges, angle);
    if (DEBUG)
        std::cout << "Angle: " << RAD2DEG(angle) << std::endl;

    // retrieve edge IDs - only 1 is expected for manifold bodies
    std::vector<int16_t> edgeIDs;
    edgeIDs.reserve(sharedEdges.Extent());
    for (int16_t i = 1; i <= sharedEdges.Extent(); ++i)
    {
        for (int16_t j = 1; j <= eMap.Extent(); ++j)
        {
            if (sharedEdges(i).IsSame(eMap(j)))
            {
                edgeIDs.push_back(j);
                break;
            }
        }
    }

    // insert edge category
    if ((FeatureAngleType_Convex == angleType ||
        FeatureAngleType_SmoothConvex == angleType) &&
        convexEdges )
    {
        InsertEdges(sharedEdges, edgeIDs, convexEdges);
    }
    else if ((FeatureAngleType_Concave == angleType ||
            FeatureAngleType_SmoothConcave == angleType) &&
            concaveEdges)
    {
        InsertEdges(sharedEdges, edgeIDs, concaveEdges);
    }
    else if (FeatureAngleType_Smooth == angleType && smoothC1Edges)
    {
        InsertEdges(sharedEdges, edgeIDs, smoothC1Edges);
    }
    else if (undefinedEdges)
    {
        InsertEdges(sharedEdges, edgeIDs, undefinedEdges);
    }

    return RAD2DEG(angle);
}


bool AS_FeatureRecognition::ClassifyEdges(std::map<int, TopoDS_Edge>* undefinedEdges)
{
    if (!this->m_AAG) return false;

    Handle(asiAlgo_AAGRandomIterator) fitr = new asiAlgo_AAGRandomIterator(this->m_AAG);
    // loop over faces
    for (; fitr->More(); fitr->Next())
    {
        auto faceID = fitr->GetFaceId();
        const TopoDS_Face& currentFace = this->m_AAG->GetFace(faceID);
        TColStd_PackedMapOfInteger neighbourIDs;
        fitr->GetNeighbors(neighbourIDs);

        // loop over neighbors
        for (TColStd_MapIteratorOfPackedMapOfInteger nfitr(neighbourIDs); nfitr.More(); nfitr.Next())
        {
            auto neighbourFaceID = nfitr.Key();
            const TopoDS_Face& neighbourFace = this->m_AAG->GetFace(neighbourFaceID);

            // calculate the angle between the two faces and accordingly classify the edge type
            auto angle = this->ComputeDihedralAngle(currentFace,
                                                    neighbourFace,
                                                    true,
                                                    1e-3,
                                                    this->m_ConvexEdges.get(),
                                                    this->m_ConcaveEdges.get(),
                                                    this->m_SmoothC1Edges.get(),
                                                    undefinedEdges );
        }
    }

    if (DEBUG)
    {
        std::cout << "\nEdge-Types:" << std::endl;
        std::cout << "Convex edges = " << this->m_ConvexEdges->size() << std::endl;
        std::cout << "Concave edges = " << this->m_ConcaveEdges->size() << std::endl;
        std::cout << "Smooth (C1) edges = " << this->m_SmoothC1Edges->size() << std::endl;
        if (undefinedEdges)
            std::cout << "Undefined edges = " << undefinedEdges->size() << std::endl;
    }
    return true;
}


bool AS_FeatureRecognition::IdentifySlots(std::vector<std::array<float, 3>>& slotDimensions)
{
    bool status = true;
    if (this->m_ConcaveEdges->empty()) return status;

    TopTools_IndexedMapOfShape vMap, eMap, fMap;
    TopTools_IndexedDataMapOfShapeListOfShape veMap, efMap;
    status &= this->m_CADPart->QueryTopologyMaps(vMap, eMap, fMap, veMap, efMap);
    if (!status) return status;

    // shared faces are the faces common to the concave edges
    // other faces are the ones not shared across the concave edges
    // shared faces and other faces are mutually exclusive
    std::cout << "\nSlot-Features:" << std::endl;
    slotDimensions.reserve(this->m_ConcaveEdges->size() + 1);
    TColStd_PackedMapOfInteger sharedFaceIDs;
    std::map<int16_t, gp_Dir> otherFaceIndexNormalMap;
    int16_t count = 0, initFace1 = -1, initFace2 = -1;
    for (const auto& eitr : *this->m_ConcaveEdges)
    {
        BRepAdaptor_Curve curve(eitr.second);
        // rectangular slots have linear concave edge
        if (GeomAbs_CurveType::GeomAbs_Line != curve.GetType())
            continue;

        // slots shared across the same face should be counted only once
        // only 2 faces per edge are counted assuming manifold geometry
        auto sharedFaces = efMap.FindFromIndex(eitr.first);
        const TopoDS_Face& face1 = TopoDS::Face(sharedFaces.First());
        const TopoDS_Face& face2 = TopoDS::Face(sharedFaces.Last());
        auto face1ID = fMap.FindIndex(face1);
        auto face2ID = fMap.FindIndex(face2);

        auto commonFaceFound = false;
        int8_t otherFace = 0;
        gp_Dir otherFaceNormal;
        auto foundFace1 = sharedFaceIDs.Contains(face1ID);
        auto foundFace2 = sharedFaceIDs.Contains(face2ID);
        if (foundFace1 || foundFace2)
        {
            if (foundFace1)
            {
                sharedFaceIDs.Add(face1ID);
                if (sharedFaceIDs.Contains(initFace2))
                {
                    sharedFaceIDs.Remove(initFace2);
                    otherFaceIndexNormalMap.erase(otherFaceIndexNormalMap.find(initFace1));
                }
            }
            else
            {
                otherFace = 1;
                otherFaceNormal = this->m_CADPart->ComputePlaneNormal(sharedFaces.First());
            }

            if (foundFace2)
            {
                sharedFaceIDs.Add(face2ID);
                if (sharedFaceIDs.Contains(initFace1))
                {
                    sharedFaceIDs.Remove(initFace1);
                    otherFaceIndexNormalMap.erase(otherFaceIndexNormalMap.find(initFace2));
                }
            }
            else
            {
                otherFace = 2;
                otherFaceNormal = this->m_CADPart->ComputePlaneNormal(sharedFaces.Last());
            }
        }
        else
        {
            initFace1 = face1ID, initFace2 = face2ID;
            gp_Dir normal1 = this->m_CADPart->ComputePlaneNormal(sharedFaces.First());
            gp_Dir normal2 = this->m_CADPart->ComputePlaneNormal(sharedFaces.Last());
            otherFaceIndexNormalMap.emplace(std::make_pair(initFace1, normal1));
            otherFaceIndexNormalMap.emplace(std::make_pair(initFace2, normal2));
            sharedFaceIDs.Add(initFace1);
            sharedFaceIDs.Add(initFace2);
        }

        if (count)
        {
            // presently common slots are identified only by faces having oppositel directed normals
            // this condition is not very robust - ideally we need to identify a common convex edge 
            // which terminates on the oppositely directed face vertices
            for (const auto& fitr : otherFaceIndexNormalMap)
            {
                if (fitr.second.IsOpposite(otherFaceNormal, Precision::Angular()))
                {
                    commonFaceFound = true;
                    break;
                }
            }
            if (1 == otherFace)
                otherFaceIndexNormalMap.emplace(std::make_pair(face1ID, otherFaceNormal));
            else if (2 == otherFace)
                otherFaceIndexNormalMap.emplace(std::make_pair(face2ID, otherFaceNormal));

            if (commonFaceFound)
                continue;
        }

        // 1 concave edge
        float lengths[3]{ 0.0f };
        int16_t edgeIDs[3]{ -1 };
        edgeIDs[0] = eitr.first;
        GProp_GProps edgeProps;
        BRepGProp::LinearProperties(eitr.second, edgeProps);
        lengths[0] = static_cast<float>(edgeProps.Mass());

        for (int16_t i = 0; i < 2; ++i)
        {
            // j = 0: original concave edge, j = 1, 2: convex edges
            int16_t j = 1;
            ShapeAnalysis_Edge edgeAnalysis;
            TopoDS_Vertex vertex;
            // choose the vertex which has 2 convex edges - if 1st vertex doesn't work, use 2nd
            if (!i)
                vertex = edgeAnalysis.FirstVertex(eitr.second);
            else
                vertex = edgeAnalysis.LastVertex(eitr.second);

            const TopTools_ListOfShape& incidentEdges = veMap.FindFromKey(vertex);
            for (TopTools_ListIteratorOfListOfShape ieitr(incidentEdges); ieitr.More(); ieitr.Next())
            {
                const TopoDS_Edge& edge = TopoDS::Edge(ieitr.Value());
                // use only convex edges here - if the current edge is concave, use the other vertex
                if (!edge.IsSame(eitr.second))
                {
                    if (this->m_ConcaveEdges->find(eMap.FindIndex(edge)) != this->m_ConcaveEdges->end())
                        break;
                }
                else
                    continue;

                BRepGProp::LinearProperties(edge, edgeProps);
                lengths[j] = static_cast<float>(edgeProps.Mass());
                edgeIDs[j] = eMap.FindIndex(edge);
                // if edge triplet is found, iterate over the next concave edge of outer loop (eitr)
                if (3 == ++j)
                    break;
            }

            if (3 == j)
                break;
        }

        std::cout << "Slot" << ++count << ":" << std::endl;
        std::cout << "Edge-ID: " << (edgeIDs[0] = eitr.first) << ", Length: " << lengths[0] << std::endl;
        for (int16_t j = 1; j < 3; ++j)
            std::cout << "Edge-ID: " << edgeIDs[j] << ", Length: " << lengths[j] << std::endl;
        slotDimensions.emplace_back(std::move(std::array<float, 3>\
            ({ lengths[0], lengths[1], lengths[2] })));
    }

    slotDimensions.shrink_to_fit();
    return status;
}


bool AS_FeatureRecognition::FindHoleFaces()
{
    asiAlgo_RecognizeDrillHoles drilledHoles(this->m_AAG);
    bool status = drilledHoles.Perform();
    if (!status) return status;
    TColStd_PackedMapOfInteger results = drilledHoles.GetResultIndices();
    if (results.IsEmpty())
    {
        drilledHoles.SetHardFeatureMode(true);
        drilledHoles.Perform();
        results = drilledHoles.GetResultIndices();
    }
    this->m_HoleFaces->Assign(results);
    return status;
}


bool AS_FeatureRecognition::IdentifyDrilledHoles(std::vector<std::pair<float, float>>& holeDimensions)
{
    // Geom_ElementarySurface is base class of Geom_CylindricalSurface, Geom_ConicalSurface
    // it used instead of gp_Cone, gp_Cylinder since they don't have a common base class
    bool status = true;
    if (this->m_HoleFaces->IsEmpty()) return status;

    std::cout << "\nHole-Features:" << std::endl;
    holeDimensions.reserve(this->m_HoleFaces->Extent());
    std::vector<std::tuple<float, gp_Ax1, gp_Pnt>> uniqueHoles;
    uniqueHoles.reserve(this->m_HoleFaces->Extent());
    for (TColStd_MapIteratorOfPackedMapOfInteger hfitr(*this->m_HoleFaces); hfitr.More(); hfitr.Next())
    {
        auto faceID = hfitr.Key();
        const TopoDS_Face& face = this->m_AAG->GetFace(faceID);
        BRepAdaptor_Surface surface(face);
        ConeCylinderDimensions quadricDims;
        bool uniqueHole = this->m_CADPart->IdentifyUniqueDegenerateQuadrics(surface, uniqueHoles);

        if (uniqueHole)
        {
            this->m_CADPart->ComputeDegenerateQuadricDimensions(surface, &quadricDims);
            auto radius = quadricDims.radius, height = quadricDims.height, angle = quadricDims.angle;
            std::cout << "Face-ID: " << faceID << ", (Radius = " << radius << ", Height = " << height;
            if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
                std::cout << ", Semi-Angle = " << angle;
            std::cout << ")" << std::endl;
        }
    }

    holeDimensions.shrink_to_fit();
    return status;
}


bool AS_FeatureRecognition::IdentifyHandles(std::vector<std::pair<float, float>>& handleDimensions)
{
    if (this->m_ConcaveEdges->empty()) return false;

    std::cout << "\nHandle-Features:" << std::endl;
    std::map<int, std::pair<float, float>> cylindricalHandlesMap;
    std::vector<std::pair<int, float>> cylindricalHandlesVector;
    std::vector<std::tuple<float, gp_Ax1, gp_Pnt>> uniqueHandles;
    uniqueHandles.reserve(this->m_ConcaveEdges->size() + 1);
    handleDimensions.reserve(this->m_ConcaveEdges->size() + 1);

    auto endCapPlanarFaceID = -1, endCylinderFaceID = -1;
    auto endCapFound = false;
    Handle(asiAlgo_AAGRandomIterator) fitr = new asiAlgo_AAGRandomIterator(this->m_AAG);
    for (; fitr->More(); fitr->Next())
    {
        auto faceID = fitr->GetFaceId();
        const TopoDS_Face& face = this->m_AAG->GetFace(faceID);
        BRepAdaptor_Surface surface(face);
        ConeCylinderDimensions quadricDims;

        if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
        {
            if (this->m_HoleFaces->Contains(faceID))
                continue;

            gp_Cylinder cylSurface = surface.Cylinder();
            bool uniqueCylinder = this->m_CADPart->IdentifyUniqueDegenerateQuadrics(surface, uniqueHandles);
            if (uniqueCylinder)
            {
                // cylinder center is on the cylinder axis, not on the cylindrical surface
                this->m_CADPart->ComputeDegenerateQuadricDimensions(surface, &quadricDims);
                auto handleRadius = quadricDims.radius, handleHeight = quadricDims.height;
                gp_Pnt handleCenter = quadricDims.center;

                std::cout << "Face-ID: " << faceID << \
                    ", (Radius = " << handleRadius << ", Height = " << handleHeight << ")" << std::endl;
                handleDimensions.emplace_back(std::make_pair(handleRadius, handleHeight));
                // center component is the cylinder center projected on the cylinder axis
                // it is a scalar and used to establish cylinder adjacencies
                float centerComponent = fabs(gp_Vec(handleCenter.XYZ()).Dot(cylSurface.Axis().Direction()));
                cylindricalHandlesMap.emplace(std::make_pair(std::move(faceID), \
                    std::make_pair(std::move(handleRadius), std::move(centerComponent))));
            }
        }
    }

    // adjacencies are established by sorting the cylindrical handles according to their center component
    std::cout << "Adjacencies:" << std::endl;
    Utility::SortMapByValue(cylindricalHandlesMap, cylindricalHandlesVector);
    for (int16_t i = 0, count = 0; i < cylindricalHandlesVector.size() - 1; ++i)
    {
        auto itr1 = cylindricalHandlesMap.find(cylindricalHandlesVector.at(i).first);
        auto itr2 = cylindricalHandlesMap.find(cylindricalHandlesVector.at(i + 1).first);
        auto diff = fabs(itr1->second.first - itr2->second.first);
        std::cout << "(" << itr1->first << ", " << itr2->first << "): " << "Differential-Radius = " << diff << std::endl;
    }

    return true;
}


bool AS_FeatureRecognition::RunFeatureRecognition(  const std::string& filePath,
                                                    int16_t iteration   )
{
    bool status = true;

    // constructor call required only when directory path is provided
    if (!this->m_CADPart)
        *this = AS_FeatureRecognition(filePath.c_str());
    if (this->Initialize())
    {
        auto fileName = Utility::FileName(filePath.c_str());
        if (iteration) std::cout << "\n\n";
        std::cout << "Part: " << fileName << std::endl;
        auto cadPart = this->GetCADPart();
        assert(!cadPart->GetShape()->IsNull());
        if (cadPart->Compute(false))
        {
            std::vector<std::array<float, 3>> slotDimensions;
            std::vector<std::pair<float, float>> holeDimensions, handleDimensions;
            auto partCategory = this->IdentifyPartCategory();
            if (AS_FeatureRecognition::PartCategory::MILLING == partCategory)
            {
                std::cout << "\nFeature Recognition:" << std::endl;
                status &= this->IdentifySlots(slotDimensions);
                if (cadPart->GetGenus() || cadPart->GetHoles())
                    status &= this->IdentifyDrilledHoles(holeDimensions);
            }
            else if (AS_FeatureRecognition::PartCategory::TURNING == partCategory)
            {
                std::cout << "\nFeature Recognition:" << std::endl;
                status &= this->IdentifyHandles(handleDimensions);
                if (cadPart->GetGenus() || cadPart->GetHoles())
                    status &= this->IdentifyDrilledHoles(holeDimensions);
            }
        }
    }
    else
        status = false;

    return status;
}
