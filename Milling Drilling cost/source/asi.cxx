#include <iomanip>
// project
#include <asi.hxx>
#include <common.hxx>
#include <occ.hxx>
// 3rd-party: ASI
#include <asiAlgo.h>
#include <asiAlgo_AAGIterator.h>
#include <asiAlgo_CheckDihedralAngle.h>
#include <asiAlgo_CheckValidity.h>
#include <asiAlgo_ConvertCanonical.h>
#include <asiAlgo_ConvertCanonicalSummary.h>
#include <asiAlgo_IGES.h>
#include <asiAlgo_InvertShells.h>
#include <asiAlgo_ReadSTEPWithMeta.h>
#include <asiAlgo_RecognizeDrillHoles.h>
#include <asiAlgo_RecognizeShafts.h>
#include <asiAlgo_STEP.h>
#include <asiAlgo_Utils.h>
#include <asiData.h>
#include <asiEngine.h>
#include <ActData.h>
// 3rd-party: OCC
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GProp_GProps.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <ShapeCustom.hxx>
#include <TopExp_Explorer.hxx>
#include <TopOpeBRepBuild_Tools.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>


AS_FeatureRecognition::AS_FeatureRecognition(   const char* cadFilePath,
                                                TopoDS_Shape* shape   )
{
    this->m_CADPart = std::make_shared<OCC_Computation>(cadFilePath, shape);
    this->m_ConvexEdges = std::make_shared<std::map<int16_t, TopoDS_Edge>>();
    this->m_ConcaveEdges = std::make_shared<std::map<int16_t, TopoDS_Edge>>();
    this->m_SmoothC1Edges = std::make_shared<std::map<int16_t, TopoDS_Edge>>();
    this->m_HoleEdges = std::make_shared<std::map<int16_t, TopoDS_Edge>>();
    this->m_Threads = std::make_shared<std::map<int16_t, std::pair<float, float>>>();
    this->m_HoleFaces = std::make_shared<TColStd_PackedMapOfInteger>();
    this->m_ShaftFaces = std::make_shared<TColStd_PackedMapOfInteger>();
    if (Utility::DEBUG_LEVEL >= 3) std::cout << __FUNCSIG__ << std::endl;
}


AS_FeatureRecognition::~AS_FeatureRecognition()
{
    if (Utility::DEBUG_LEVEL >= 3) std::cout << __FUNCSIG__ << std::endl;
}


void AS_FeatureRecognition::MoveSemantics(AS_FeatureRecognition& rhs) noexcept
{
    this->m_CADPart = std::move(rhs.m_CADPart);
    this->m_ConvexEdges = std::move(rhs.m_ConvexEdges);
    this->m_ConcaveEdges = std::move(rhs.m_ConcaveEdges);
    this->m_SmoothC1Edges = std::move(rhs.m_SmoothC1Edges);
    this->m_HoleEdges = std::move(rhs.m_HoleEdges);
    this->m_Threads = std::move(rhs.m_Threads);
    this->m_HoleFaces = std::move(rhs.m_HoleFaces);
    this->m_ShaftFaces = std::move(rhs.m_ShaftFaces);
    this->m_AAG = std::move(rhs.m_AAG);
}


AS_FeatureRecognition::AS_FeatureRecognition(AS_FeatureRecognition&& rhs) noexcept
{
    this->MoveSemantics(rhs);
}


AS_FeatureRecognition& AS_FeatureRecognition::operator=(AS_FeatureRecognition&& rhs) noexcept
{
    // check memory address only, not the pointer contents
    if (this != &rhs)
        this->MoveSemantics(rhs);
    return *this;
}


bool AS_FeatureRecognition::STEPReader()
{
    // healing of step files is set to true
    Handle(asiAlgo_STEP) stp = new asiAlgo_STEP;
    return stp->Read(TCollection_AsciiString(this->m_CADPart->GetFilePath().c_str()), true, *this->m_CADPart->GetShape().get());
    //Handle(asiAlgo_ReadSTEPWithMeta) reader = new asiAlgo_ReadSTEPWithMeta;
    //bool status = reader->Perform(this->m_FilePath.c_str());
}


bool AS_FeatureRecognition::IGESReader()
{
    Handle(asiAlgo_IGES) igs = new asiAlgo_IGES(ActAPI_ProgressEntry());
    bool success = igs->Read(TCollection_AsciiString(this->m_CADPart->GetFilePath().c_str()), *this->m_CADPart->GetShape().get());

    if (success)
    {
        // topo-geom maps are needed for computation of sewing tolerance (ConvertIGESToSolid)
        if (!this->m_CADPart->GenerateTopologyGeometryMaps()) return false;
        success &= this->m_CADPart->ConvertIGESToSolid();
        if (!success) return success;
        success &= this->InvertShellOrientation();
        if (!success) return success;
        // updation of topo-geom maps is needed after converting surfaces to solid
        // and shell inversion (if executed)
        // in the surface form, there are no adjacencies established (veMap, efMap)
        if (!this->m_CADPart->UpdateTopologyGeometryMaps()) return false;
        // canonical conversion works only on interpolated geometry (Bezier/BSpline) and not
        // on procedural geometry, hence this is an intermediate step
        std::pair<int16_t, int16_t> convertedEntities;
        TopoDS_Shape modifiedShape;
        if (Utility::DEBUG_LEVEL >= 2)
        {
            this->m_CADPart->ComputeCurveSurfaceGeometryMaps();
            convertedEntities = this->m_CADPart->CountProceduralGeometryEntities();
        }
#if 1
        // handles both curves & surfaces
        modifiedShape = ShapeCustom::ConvertToBSpline(*this->m_CADPart->GetShape(), true, true, true);
#else
        auto convertedEntities = this->m_CADPart->ConvertProceduralGeometryToBSplineGeometry();
#endif
        this->m_CADPart->SetShape(std::make_shared<TopoDS_Shape>(modifiedShape));
        if (!this->m_CADPart->UpdateTopologyGeometryMaps()) return false;
        if (Utility::DEBUG_LEVEL >= 2)
        {
            this->m_CADPart->ComputeCurveSurfaceGeometryMaps();
            auto temp = this->m_CADPart->CountProceduralGeometryEntities();
            convertedEntities.first -= temp.first;
            convertedEntities.second -= temp.second;
            if (convertedEntities.first > 0 || convertedEntities.second > 0)
            {
                std::cout << "Procedural Geometry Conversions: Curves: " << convertedEntities.first << \
                    ", Surfaces: " << convertedEntities.second << std::endl;
            }
        }
        // conversion to canonical geometry is required for fixing redundant topology and also for
        // downstream feature recognition
        success &= this->ConvertToCanonicalGeometry();
        if (!success) return success;
        if (!this->m_CADPart->UpdateTopologyGeometryMaps()) return false;
        // updation of topo-geom maps is needed after converting freeform geometry to canonical geometry
        // (analytical form)
        auto reducedEntities = this->m_CADPart->FixRedundantTopology();
        if (Utility::NonZeroContainerCount(reducedEntities) > 0)
        {
            if (Utility::DEBUG_LEVEL >= 2)
            {
                std::cout << "Topological Reductions: Vertices: " << reducedEntities[0] << \
                    ", Edges: " << reducedEntities[1] << ", Faces: " << reducedEntities[2] << "\n";
            }
            if (!this->m_CADPart->UpdateTopologyGeometryMaps()) return false;
#if 1
            // TBD : after fixing redundant topology, the altered edges get converted back to interpolated geometry
            // ideally conversion to canonical geometry should not be required
            success &= this->ConvertToCanonicalGeometry();
            if (!success) return success;
            if (!this->m_CADPart->UpdateTopologyGeometryMaps()) return false;
#endif
        }
    }

    return success;
}


bool AS_FeatureRecognition::GenerateAAG()
{
    bool success = true;
    const auto defaultPrecision = std::cout.precision();
    if (this->m_CADPart->GetShape()->IsNull())
    {
        auto fileExtn = Utility::FileExtn(this->m_CADPart->GetFilePath().c_str());
        if ("step" == fileExtn || "stp" == fileExtn)
            success &= this->STEPReader();
        else if ("iges" == fileExtn || "igs" == fileExtn)
            success &= this->IGESReader();
    }
    if (!this->m_CADPart->GetTopoGeomMapsGenerated())
        success &= this->m_CADPart->GenerateTopologyGeometryMaps();

    if (success)
    {
        // attribute adjacency graph
        this->m_AAG = new asiAlgo_AAG(*this->m_CADPart->GetShape());
        if (defaultPrecision != std::cout.precision())
        {
            std::cout << std::setprecision(defaultPrecision);
            std::cout << std::defaultfloat;
        }
#if 0
        //this->m_AAG->Dump(std::cout);
        //this->m_CADPart->WriteFile();                       // for testing only
        this->m_CADPart->ComputePartGeometry();             // for testing only
#endif
    }

    return success;
}


bool AS_FeatureRecognition::Initialize()
{
    bool status = this->GenerateAAG();
    if (!status) return status;
    status &= this->ClassifyEdges();
    status &= this->ComputeHoleFaces();
    status &= this->ComputeShaftFaces();
    //...
    return status;
}


bool AS_FeatureRecognition::InvertShellOrientation()
{
    bool status = true;
    // test for negative volume and if it tests positive, invert the shells
    auto topoShape = this->m_CADPart->GetShape();
    GProp_GProps massProps;
    BRepGProp::VolumeProperties(*topoShape, massProps);
    if (massProps.Mass() < Precision::Confusion())
    {
        asiAlgo_InvertShells inverter(*topoShape);
        status &= inverter.Perform();
        if (status)
            this->m_CADPart->SetShape(std::make_shared<TopoDS_Shape>(inverter.GetResult()));
    }
    return status;
}


bool AS_FeatureRecognition::ConvertToCanonicalGeometry()
{
    auto cadPart = this->m_CADPart;
    if (cadPart->GetShape()->IsNull()) return false;

    // stitching tolerance (OCC_Computation::m_Tolerance) is generally on the higher end,
    // hence shape tolerance of combined vertex, edge, face is used
    ShapeAnalysis_ShapeTolerance shapeTolerance;
    auto maxTol = shapeTolerance.Tolerance(*cadPart->GetShape(), 1);
    asiAlgo_ConvertCanonical converter;
    TopoDS_Shape modifiedShape = converter.Perform(*cadPart->GetShape(), maxTol);
    asiAlgo_CheckValidity checker;
    bool status = checker.CheckBasic(modifiedShape);
    if (!status) return status;

    // feature recognition section of the code is not supposed to modify the original CAD shape
    // but in case of IGES files, feature recognition is easier with canonical (analytical) shapes
    this->m_CADPart->SetShape(std::make_shared<TopoDS_Shape>(modifiedShape));
    asiAlgo_ConvertCanonicalSummary summary = converter.GetSummary();
    summary.isValid = true;
    if (Utility::DEBUG_LEVEL >= 2)
    {
        std::cout << "Canonical Geometry Conversions:";
        asiAlgo_ConvertCanonicalSummary::ToJSON(summary, 0, std::cout);
        std::cout << "\n";
    }
    return status;
}


AS_FeatureRecognition::PartCategory AS_FeatureRecognition::IdentifyPartCategory()
{
    AS_FeatureRecognition::PartCategory enumObj = AS_FeatureRecognition::PartCategory::UNDEFINED;

    std::map<int16_t, GeomAbs_CurveType> etMap;
    std::map<int16_t, GeomAbs_SurfaceType> ftMap;
    bool status = this->m_CADPart->QueryGeometryMaps(etMap, ftMap);
    if (!status) return enumObj;
    status &= this->ComputeHoleEdges();
    if (!status) return enumObj;

    // this criteria will change as more CAD part categories are introduced
    for (const auto& itr : *this->m_ConcaveEdges)
    {
        // find only those edges which are not part of holes
        if (this->m_HoleEdges->end() != this->m_HoleEdges->find(itr.first))
            continue;
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
                                                    std::map<int16_t, TopoDS_Edge>* convexEdges,
                                                    std::map<int16_t, TopoDS_Edge>* concaveEdges,
                                                    std::map<int16_t, TopoDS_Edge>* smoothC1Edges,
                                                    std::map<int16_t, TopoDS_Edge>* undefinedEdges  ) const
{
    auto InsertEdges = [](  const TopTools_IndexedMapOfShape& sharedEdges,
                            const std::vector<int16_t>& edgeIDs,
                            std::map<int16_t, TopoDS_Edge>* edgeType) -> void
    {
        for (int16_t idx = 1; idx <= sharedEdges.Extent(); ++idx)
            edgeType->insert(std::make_pair(edgeIDs.at(static_cast<int64_t>(idx) - 1), TopoDS::Edge(sharedEdges.FindKey(idx))));
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
#if 0
    if (Utility::DEBUG_LEVEL >= 3)
        std::cout << "Angle: " << RAD2DEG(angle) << std::endl;
#endif

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


bool AS_FeatureRecognition::ClassifyEdges(std::map<int16_t, TopoDS_Edge>* undefinedEdges)
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

    if (Utility::DEBUG_LEVEL >= 3)
    {
        std::cout << "\nEdge-Types:" << std::endl;
        std::cout << "Convex edges = " << this->m_ConvexEdges->size() << std::endl;
        std::cout << "Concave edges = " << this->m_ConcaveEdges->size() << std::endl;
        std::cout << "Smooth (C1) edges = " << this->m_SmoothC1Edges->size() << std::endl;
        if (undefinedEdges)
            std::cout << "Undefined edges = " << undefinedEdges->size() << std::endl;
        std::cout << "Concave-Edge-IDs:\n";
        for (const auto& itr : *this->m_ConcaveEdges)
        {
            std::cout << itr.first << ": ";
            std::cout << OCC_Computation::CurveGeometry(BRepAdaptor_Curve(itr.second).GetType()) << std::endl;
        }
        std::cout << "\n";
    }
    return true;
}


bool AS_FeatureRecognition::IdentifySlots(std::vector<std::array<float, 3>>& slotDimensions) const
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
    if (Utility::DEBUG_LEVEL >= 1)
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
                otherFaceNormal = OCC_Computation::ComputePlaneNormal(sharedFaces.First());
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
                otherFaceNormal = OCC_Computation::ComputePlaneNormal(sharedFaces.Last());
            }
        }
        else
        {
            initFace1 = face1ID, initFace2 = face2ID;
            gp_Dir normal1 = OCC_Computation::ComputePlaneNormal(sharedFaces.First());
            gp_Dir normal2 = OCC_Computation::ComputePlaneNormal(sharedFaces.Last());
            otherFaceIndexNormalMap.emplace(initFace1, normal1); // std::make_pair()
            otherFaceIndexNormalMap.emplace(initFace2, normal2); // std::make_pair()
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
                otherFaceIndexNormalMap.emplace(face1ID, otherFaceNormal); // std::make_pair()
            else if (2 == otherFace)
                otherFaceIndexNormalMap.emplace(face2ID, otherFaceNormal); // std::make_pair()

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

        slotDimensions.emplace_back(std::move(std::array<float, 3>\
            ({ lengths[0], lengths[1], lengths[2] })));
        if (Utility::DEBUG_LEVEL >= 1)
        {
            std::cout << "Slot" << ++count << ":" << std::endl;
            std::cout << "Edge-ID: " << (edgeIDs[0] = eitr.first) << ", Length: " << lengths[0] << std::endl;
            for (int16_t j = 1; j < 3; ++j)
                std::cout << "Edge-ID: " << edgeIDs[j] << ", Length: " << lengths[j] << std::endl;
        }
    }

    slotDimensions.shrink_to_fit();
    return status;
}


void AS_FeatureRecognition::PrintFaceIDGeometry(bool featureCategory) const
{
    std::shared_ptr<TColStd_PackedMapOfInteger> faces;
    if (featureCategory)
        faces = this->m_HoleFaces;
    else
        faces = this->m_ShaftFaces;
    for (TColStd_MapIteratorOfPackedMapOfInteger fitr(*faces); fitr.More(); fitr.Next())
    {
        auto faceID = fitr.Key();
        std::cout << faceID << ": ";
        std::cout << OCC_Computation::SurfaceGeometry(BRepAdaptor_Surface(this->m_AAG->GetFace(faceID)).GetType()) \
            << std::endl;
    }
}


bool AS_FeatureRecognition::ComputeHoleFaces()
{
    asiAlgo_RecognizeDrillHoles drilledHoles(this->m_AAG, true);
    bool status = drilledHoles.Perform();
    if (!status) return status;
    TColStd_PackedMapOfInteger results = drilledHoles.GetResultIndices();
    this->m_HoleFaces->Assign(results);
    if (Utility::DEBUG_LEVEL >= 3)
    {
        std::cout << "Hole-Faces = " << this->m_HoleFaces->Extent() << std::endl;
        std::cout << "Hole-Face-IDs:\n";
        this->PrintFaceIDGeometry(true);
        std::cout << "\n";
    }
    return status;
}


bool AS_FeatureRecognition::ComputeHoleEdges()
{
    bool status = true;
    if (this->m_HoleFaces->IsEmpty()) return status;

    TopTools_IndexedMapOfShape vMap, eMap, fMap;
    TopTools_IndexedDataMapOfShapeListOfShape veMap, efMap;
    status &= this->m_CADPart->QueryTopologyMaps(vMap, eMap, fMap, veMap, efMap);
    if (!status) return status;

    // every cylindrical/conical face should have a helical edge for threaded holes
    // if some faces do not have helical edge present, most likely it is not a threaded hole
    std::vector<bool> helicalEdges(this->m_HoleFaces->Extent(), false);
    // std::vector<bool>::iterator doesn't behave as expected, hence integer counter is used instead
    auto count = 0;

    std::vector<std::vector<OCC_Computation::HelixDimensions>> sortedHelices;
    std::vector<OCC_Computation::HelixDimensions> helixDimensions;
    OCC_Computation::HelixDimensions dimensions;
    helixDimensions.reserve(static_cast<int64_t>(this->m_HoleFaces->Extent()) * 5);
    for (TColStd_MapIteratorOfPackedMapOfInteger hfitr(*this->m_HoleFaces); hfitr.More(); hfitr.Next())
    {
        auto faceID = hfitr.Key();
        const TopoDS_Face& face = this->m_AAG->GetFace(faceID);
        BRepAdaptor_Surface surface(face);
        for (TopExp_Explorer exp(face, TopAbs_EDGE); exp.More(); exp.Next())
        {
            const TopoDS_Edge& edge = TopoDS::Edge(exp.Current());
            auto edgeID = eMap.FindIndex(edge);
            this->m_HoleEdges->emplace(edgeID, edge); // std::make_pair()
            // if helical edge is detected, no need to compute again for the same face
            if (count < helicalEdges.size() && !helicalEdges[count])
            {
                // helix wraps itself around a cone or cylinder
                if (GeomAbs_Cylinder == surface.GetType() ||
                    GeomAbs_Cone == surface.GetType())
                {
                    // OCC doesn't have in-built helix curve, it is categorized as BSpline or Bezier
                    BRepAdaptor_Curve curve(edge);
                    if (GeomAbs_BezierCurve == curve.GetType() ||
                        GeomAbs_BSplineCurve == curve.GetType())
                    {
                        if (OCC_Computation::IsInterpolatedCurveHelix(surface, curve, face, edge, faceID, edgeID, &dimensions))
                        {
                            // helixDimensions holds the data of all helix curves from multiple holes in the CAD model
                            // this needs to be sorted later according to their axes to compute mode of helix radius & pitch
                            helicalEdges[count] = true;
                            helixDimensions.emplace_back(std::move(dimensions));
                        }
                        else
                            helicalEdges.erase(helicalEdges.begin() + count);
                    }
                }
                else
                    helicalEdges.erase(helicalEdges.begin() + count);
            }
        }

        // increment the iterator only if the element is not deleted
        if (count < helicalEdges.size() && helicalEdges[count]) ++count;
    }

    // at least 1 interpolated curve should be classified as helix for threads to be present
    count = std::count_if(helicalEdges.begin(), helicalEdges.end(), [](auto idx) { return true == idx; });
    if (count)
    {
        this->m_ThreadsPresent = true;
        OCC_Computation::SortHelicesByAxes(helixDimensions, sortedHelices);
        for (const auto& itr : sortedHelices)
        {
            OCC_Computation::ComputeHelixModeDimensions(itr, &dimensions);
            this->m_Threads->emplace(dimensions.parentFace,\
                std::make_pair(dimensions.radius, dimensions.pitch));
        }
    }

    return status;
}


bool AS_FeatureRecognition::ComputeShaftFaces()
{
    if (!this->m_HoleFaces.get()) return false;
    asiAlgo_RecognizeShafts shafts(this->m_AAG);
    shafts.SetExcludedFaces(*this->m_HoleFaces);
    bool status = shafts.Perform();
    if (!status) return status;
    TColStd_PackedMapOfInteger results = shafts.GetResultIndices();
    if (!results.IsEmpty())
        this->m_ShaftFaces->Assign(results);

    if (Utility::DEBUG_LEVEL >= 3 && !this->m_ShaftFaces->IsEmpty())
    {
        std::cout << "Shaft-Faces = " << this->m_ShaftFaces->Extent() << std::endl;
        std::cout << "Shaft-Face-IDs:\n";
        this->PrintFaceIDGeometry(false);
        std::cout << "\n";
    }
    return status;
}


bool AS_FeatureRecognition::AdjacentHelicoidPresent(t_topoId faceID) const
{
    bool threadedFacePresent = false;
    if (this->m_ThreadsPresent)
    {
        TColStd_PackedMapOfInteger neighbourIDs;
        neighbourIDs = this->m_AAG->GetNeighborsThruVerts(faceID);
        for (TColStd_MapIteratorOfPackedMapOfInteger nfitr(neighbourIDs); nfitr.More(); nfitr.Next())
        {
            // check for neighbouring face
            auto neighbourFaceID = nfitr.Key();
            const TopoDS_Face& neighbourFace = this->m_AAG->GetFace(neighbourFaceID);
            BRepAdaptor_Surface neighbourSurface(neighbourFace);
            if (GeomAbs_BSplineSurface == neighbourSurface.GetType() ||
                GeomAbs_BezierSurface == neighbourSurface.GetType())
            {
                threadedFacePresent = true;
                break;
            }
        }

        if (!threadedFacePresent)
        {
            for (TColStd_MapIteratorOfPackedMapOfInteger nfitr(neighbourIDs); nfitr.More(); nfitr.Next())
            {
                auto neighbourFaceID = nfitr.Key();
                const TopoDS_Face& neighbourFace = this->m_AAG->GetFace(neighbourFaceID);
                for (TopExp_Explorer exp(neighbourFace, TopAbs_EDGE); exp.More(); exp.Next())
                {
                    // check for edges of neighbouring face
                    const TopoDS_Edge& neighbourEdge = TopoDS::Edge(exp.Current());
                    BRepAdaptor_Curve neighbourCurve(neighbourEdge);
                    if (GeomAbs_BSplineCurve == neighbourCurve.GetType() ||
                        GeomAbs_BezierCurve == neighbourCurve.GetType())
                    {
                        threadedFacePresent = true;
                        break;
                    }
                }
                if (threadedFacePresent) break;
            }
        }
    }
    return threadedFacePresent;
}


bool AS_FeatureRecognition::IdentifyHoles(  bool drilledHoles,
                                            std::vector<std::pair<float, float>>& holeDimensions) const
{
    // Geom_ElementarySurface is base class of Geom_CylindricalSurface, Geom_ConicalSurface
    // it used instead of gp_Cone, gp_Cylinder since they don't have a common base class
    bool status = true;
    if (this->m_HoleFaces->IsEmpty()) return status;

    if (Utility::DEBUG_LEVEL >= 1)
    {
        if (drilledHoles)
            std::cout << "\nDrilled-Hole-Features:" << std::endl;
        else
            std::cout << "\nThreaded-Hole-Features:" << std::endl;
    }
    holeDimensions.reserve(this->m_HoleFaces->Extent());
    std::vector<std::tuple<float, gp_Ax1, gp_Pnt>> uniqueHoles;
    uniqueHoles.reserve(this->m_HoleFaces->Extent());
    for (TColStd_MapIteratorOfPackedMapOfInteger hfitr(*this->m_HoleFaces); hfitr.More(); hfitr.Next())
    {
        auto faceID = hfitr.Key();
        if (drilledHoles)
        {
            if (this->AdjacentHelicoidPresent(faceID))
                continue;

            const TopoDS_Face& face = this->m_AAG->GetFace(faceID);
            BRepAdaptor_Surface surface(face);
            OCC_Computation::ConeCylinderDimensions quadricDims;
            bool uniqueHole = OCC_Computation::IdentifyUniqueDegenerateQuadrics(surface, uniqueHoles);
            if (uniqueHole)
            {
                OCC_Computation::ComputeDegenerateQuadricDimensions(surface, &quadricDims);
                quadricDims.center = std::get<2>(uniqueHoles.back());
                auto radius = quadricDims.radius, height = quadricDims.height, angle = quadricDims.angle;
                holeDimensions.emplace_back(radius, height); // std::make_pair()
                if (Utility::DEBUG_LEVEL >= 1)
                {
                    std::cout << "Face-ID: " << faceID << ", (Radius = " << radius << ", Height = " << height;
                    if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
                        std::cout << ", Semi-Angle = " << angle;
                    std::cout << ")" << std::endl;
                }
            }
        }
        else
        {
            auto tfitr = this->m_Threads->find(faceID);
            if (this->m_Threads->end() != tfitr)
            {
                auto radius = tfitr->second.first;
                auto pitch = tfitr->second.second;
                holeDimensions.emplace_back(radius, pitch);
                if (Utility::DEBUG_LEVEL >= 1)
                    std::cout << "Face-ID: " << faceID << ", (Radius = " << \
                    radius << ", Pitch = " << pitch << ")" << std::endl;
            }
        }
    }

    holeDimensions.shrink_to_fit();
    return status;
}


bool AS_FeatureRecognition::IdentifyHandles(std::vector<std::pair<float, float>>& handleDimensions) const
{
    if (this->m_ConcaveEdges->empty()) return false;

    if (Utility::DEBUG_LEVEL >= 1)
        std::cout << "\nHandle-Features:" << std::endl;
    handleDimensions.reserve(this->m_ConcaveEdges->size() + 1);
    std::map<int16_t, std::pair<float, float>> cylindricalHandlesMap;
    std::vector<std::pair<int16_t, float>> cylindricalHandlesVector;
    std::vector<std::tuple<float, gp_Ax1, gp_Pnt>> uniqueHandles;
    uniqueHandles.reserve(this->m_ConcaveEdges->size() + 1);
    for (TColStd_MapIteratorOfPackedMapOfInteger sfitr(*this->m_ShaftFaces); sfitr.More(); sfitr.Next())
    {
        auto faceID = sfitr.Key();
        const TopoDS_Face& face = this->m_AAG->GetFace(faceID);
        BRepAdaptor_Surface surface(face);
        OCC_Computation::ConeCylinderDimensions quadricDims;

        if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
        {
            if (this->AdjacentHelicoidPresent(faceID))
                continue;

            gp_Cylinder cylSurface = surface.Cylinder();
            bool uniqueShaft = OCC_Computation::IdentifyUniqueDegenerateQuadrics(surface, uniqueHandles);
            if (uniqueShaft)
            {
                // cylinder center is on the cylinder axis, not on the cylindrical surface
                OCC_Computation::ComputeDegenerateQuadricDimensions(surface, &quadricDims);
                quadricDims.center = std::get<2>(uniqueHandles.back());
                auto handleRadius = quadricDims.radius, handleHeight = quadricDims.height;
                gp_Pnt handleCenter = quadricDims.center;

                handleDimensions.emplace_back(handleRadius, handleHeight); // std::make_pair()
                // center component is the cylinder center projected on the cylinder axis
                // it is a scalar and used to establish cylinder adjacencies
                float centerComponent = fabs(gp_Vec(handleCenter.XYZ()).Dot(cylSurface.Axis().Direction()));
                cylindricalHandlesMap.emplace(std::make_pair(faceID, \
                    std::make_pair(handleRadius, centerComponent)));
                if (Utility::DEBUG_LEVEL >= 1)
                {
                    std::cout << "Face-ID: " << faceID << \
                        ", (Radius = " << handleRadius << ", Height = " << handleHeight << ")" << std::endl;
                }
            }
        }
    }

    // adjacencies are established by sorting the cylindrical handles according to their center component
    if (handleDimensions.size() >= 2)
    {
        if (Utility::DEBUG_LEVEL >= 1)
            std::cout << "Adjacencies:" << std::endl;
        Utility::SortMapByValue(cylindricalHandlesMap, cylindricalHandlesVector);
        for (int16_t i = 0, count = 0; i < cylindricalHandlesVector.size() - 1; ++i)
        {
            auto itr1 = cylindricalHandlesMap.find(cylindricalHandlesVector.at(i).first);
            auto itr2 = cylindricalHandlesMap.find(cylindricalHandlesVector.at(static_cast<int64_t>(i) + 1).first);
            auto diff = fabs(itr1->second.first - itr2->second.first);
            if (diff > 1e-3 && Utility::DEBUG_LEVEL >= 1)
                std::cout << "(" << itr1->first << ", " << itr2->first << "): " << "Differential-Radius = " << diff << std::endl;
        }
    }

    return true;
}


bool AS_FeatureRecognition::RunFeatureRecognition(  const std::string& filePath,
                                                    std::vector<std::array<float, 3>>& slotDimensions,
                                                    std::vector<std::pair<float, float>>& drilledHoleDimensions,
                                                    std::vector<std::pair<float, float>>& threadedHoleDimensions,
                                                    std::vector<std::pair<float, float>>& handleDimensions,
                                                    std::array<float, 3>& boxDimensions )
{
    bool status = true;

    // (move) assignment operator & parameterized contructor call required only when directory path is provided
    if (!this->m_CADPart)
        *this = AS_FeatureRecognition(filePath.c_str());
    auto fileName = Utility::FileName(filePath.c_str());
    if (AS_FeatureRecognition::count++) std::cout << "\n\n";
    std::cout << "Part: " << fileName << std::endl;
    if (this->Initialize())
    {
        auto cadPart = this->GetCADPart();
        assert(!cadPart->GetShape()->IsNull());
        gp_XYZ boundingBox;
        if (cadPart->Compute(false, &boundingBox))
        {
            for (int8_t i = 0; i < 3; ++i)
                boxDimensions[i] = boundingBox.Coord(i + 1);
            auto partCategory = this->IdentifyPartCategory();
            if (AS_FeatureRecognition::PartCategory::MILLING == partCategory)
            {
                if (Utility::DEBUG_LEVEL >= 1)
                    std::cout << "\nFeature Recognition:" << std::endl;
                status &= this->IdentifySlots(slotDimensions);
                if (cadPart->GetGenus() || cadPart->GetHoles())
                {
                    status &= this->IdentifyHoles(true, drilledHoleDimensions);
                    if (this->m_ThreadsPresent)
                        status &= this->IdentifyHoles(false, threadedHoleDimensions);
                }
            }
            else if (AS_FeatureRecognition::PartCategory::TURNING == partCategory)
            {
                if (Utility::DEBUG_LEVEL >= 1)
                    std::cout << "\nFeature Recognition:" << std::endl;
                status &= this->IdentifyHandles(handleDimensions);
                if (cadPart->GetGenus() || cadPart->GetHoles())
                {
                    status &= this->IdentifyHoles(true, drilledHoleDimensions);
                    if (this->m_ThreadsPresent)
                        status &= this->IdentifyHoles(false, threadedHoleDimensions);
                }
            }
            else
                status = false;
        }
    }
    else
        status = false;

    return status;
}
