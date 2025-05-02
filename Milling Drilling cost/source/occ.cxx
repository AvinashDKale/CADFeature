#include <algorithm>
#include <numeric>
// project
#include <common.hxx>
#include <occ.hxx>
// 3rd-party: OCC
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Curve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomConvert.hxx>
#include <GeomLProp_CLProps.hxx>
#include <GProp_GProps.hxx>
#include <IGESControl_Reader.hxx> 
#include <IGESControl_Writer.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <ShapeAnalysis_ShapeContents.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeCustom.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>
#include <Standard_Version.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopOpeBRepTool_ShapeExplorer.hxx>
#include <TopTools_ShapeSet.hxx>
#include <Transfer_TransientProcess.hxx>
#include <XSControl_Reader.hxx>
#include <XSControl_Writer.hxx>
#include <XSControl_TransferReader.hxx>
#include <XSControl_WorkSession.hxx>


OCC_Computation::OCC_Computation(   const char* cadFilePath,
                                    TopoDS_Shape* shape )
{
    if (cadFilePath)
        this->m_FilePath.assign(cadFilePath);
    if (shape)
        this->m_Shape = std::make_shared<TopoDS_Shape>(*shape);
    else
        this->m_Shape = std::make_shared<TopoDS_Shape>();
    this->m_AxisAlignedBoundingBox = std::make_shared<Bnd_Box>();
    this->m_OrientedBoundingBox = std::make_shared<Bnd_OBB>();
    this->m_VertexMap = std::make_shared<TopTools_IndexedMapOfShape>();
    this->m_EdgeMap = std::make_shared<TopTools_IndexedMapOfShape>();
    this->m_FaceMap = std::make_shared<TopTools_IndexedMapOfShape>();
    this->m_VertexEdgeMap = std::make_shared<TopTools_IndexedDataMapOfShapeListOfShape>();
    this->m_EdgeFaceMap = std::make_shared<TopTools_IndexedDataMapOfShapeListOfShape>();
    this->m_EdgeTypeMap = std::make_shared<std::map<int16_t, GeomAbs_CurveType>>();
    this->m_FaceTypeMap = std::make_shared<std::map<int16_t, GeomAbs_SurfaceType>>();
    if (Utility::DEBUG_LEVEL >= 2)
    {
        this->m_CurveGeometryMap = std::make_shared<std::map<GeomAbs_CurveType, int16_t>>();
        this->m_SurfaceGeometryMap = std::make_shared<std::map<GeomAbs_SurfaceType, int16_t>>();
    }
    if (Utility::DEBUG_LEVEL >= 3) std::cout << __FUNCSIG__ << std::endl;
}


OCC_Computation::~OCC_Computation()
{
    if (Utility::DEBUG_LEVEL >= 3) std::cout << __FUNCSIG__ << std::endl;
}


void OCC_Computation::MoveSemantics(OCC_Computation& rhs) noexcept
{
    this->m_FilePath = std::move(rhs.m_FilePath);
    this->m_Shape = std::move(rhs.m_Shape);
    this->m_OrientedBoundingBox = std::move(rhs.m_OrientedBoundingBox);
    this->m_AxisAlignedBoundingBox = std::move(rhs.m_AxisAlignedBoundingBox);
    this->m_VertexMap = std::move(rhs.m_VertexMap);
    this->m_EdgeMap = std::move(rhs.m_EdgeMap);
    this->m_FaceMap = std::move(rhs.m_FaceMap);
    this->m_VertexEdgeMap = std::move(rhs.m_VertexEdgeMap);
    this->m_EdgeFaceMap = std::move(rhs.m_EdgeFaceMap);
    this->m_EdgeTypeMap = std::move(rhs.m_EdgeTypeMap);
    this->m_FaceTypeMap = std::move(rhs.m_FaceTypeMap);
}


OCC_Computation::OCC_Computation(OCC_Computation&& rhs) noexcept
{
    this->MoveSemantics(rhs);
}


OCC_Computation& OCC_Computation::operator=(OCC_Computation&& rhs) noexcept
{
    // check memory address only, not the pointer contents
    if (this != &rhs)
        this->MoveSemantics(rhs);
    return *this;
}


bool OCC_Computation::GenerateTopologyGeometryMaps()
{
    if (this->m_Shape->IsNull()) return false;
    if (this->m_TopoGeomMapsGenerated) return true;

    TopExp::MapShapes(*this->m_Shape, TopAbs_VERTEX, *this->m_VertexMap);
    TopExp::MapShapes(*this->m_Shape, TopAbs_EDGE, *this->m_EdgeMap);
    TopExp::MapShapes(*this->m_Shape, TopAbs_FACE, *this->m_FaceMap);
    TopExp::MapShapesAndUniqueAncestors(*this->m_Shape, TopAbs_VERTEX, TopAbs_EDGE, *this->m_VertexEdgeMap);
    TopExp::MapShapesAndUniqueAncestors(*this->m_Shape, TopAbs_EDGE, TopAbs_FACE, *this->m_EdgeFaceMap);

    for (int16_t idx = 1; idx <= this->m_EdgeMap->Extent(); ++idx)
    {
        //BRepAdaptor_Curve curve(TopoDS::Edge((*this->m_EdgeMap)(idx)));
        this->m_EdgeTypeMap->emplace(std::make_pair(idx, BRepAdaptor_Curve(TopoDS::Edge((*this->m_EdgeMap)(idx))).GetType()));
    }
    for (int16_t idx = 1; idx <= this->m_FaceMap->Extent(); ++idx)
    {
        //BRepAdaptor_Surface surface(TopoDS::Face(faceID));
        this->m_FaceTypeMap->emplace(std::make_pair(idx, BRepAdaptor_Surface(TopoDS::Face((*this->m_FaceMap)(idx))).GetType()));
    }

    this->m_TopoGeomMapsGenerated = true;
    return true;
}


void OCC_Computation::ResetTopologyGeometryMaps()
{
    if (!this->m_TopoGeomMapsGenerated) return;
    this->m_VertexMap->Clear();
    this->m_EdgeMap->Clear();
    this->m_FaceMap->Clear();
    this->m_VertexEdgeMap->Clear();
    this->m_EdgeFaceMap->Clear();
    this->m_EdgeTypeMap->clear();
    this->m_FaceTypeMap->clear();
    this->m_TopoGeomMapsGenerated = false;
}


bool OCC_Computation::UpdateTopologyGeometryMaps()
{
    this->ResetTopologyGeometryMaps();
    return this->GenerateTopologyGeometryMaps();
}


bool OCC_Computation::QueryTopologyMaps(TopTools_IndexedMapOfShape& vertexMap,
                                        TopTools_IndexedMapOfShape& edgeMap,
                                        TopTools_IndexedMapOfShape& faceMap,
                                        TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap,
                                        TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap  ) const
{
    if (!this->m_VertexMap || this->m_VertexMap->IsEmpty()) return false;
    vertexMap = *this->m_VertexMap;
    edgeMap = *this->m_EdgeMap;
    faceMap = *this->m_FaceMap;
    vertexEdgeMap = *this->m_VertexEdgeMap;
    edgeFaceMap = *this->m_EdgeFaceMap;
    return true;
}


bool OCC_Computation::QueryGeometryMaps(std::map<int16_t, GeomAbs_CurveType>& edgeTypeMap,
                                        std::map<int16_t, GeomAbs_SurfaceType>& faceTypeMap) const
{
    if (!this->m_VertexMap || this->m_VertexMap->IsEmpty()) return false;
    edgeTypeMap = *this->m_EdgeTypeMap;
    faceTypeMap = *this->m_FaceTypeMap;
    return true;
}


bool OCC_Computation::Compute(  bool usingOCC,
                                gp_XYZ* boundingBox )
{
    bool status = true;
    if (!this->m_Shape->IsNull())
    {
        if (Utility::DEBUG_LEVEL >= 1)
            std::cout << "\nCAD Computation:\n\n";
        //if (this->m_AxisAlignedBoundingBox->SquareExtent() < Precision::SquareConfusion())
        gp_XYZ box = this->ComputeBoundingBox();
        if (boundingBox) *boundingBox = box;
        if (usingOCC && !this->m_TopoGeomMapsGenerated)
            status &= this->GenerateTopologyGeometryMaps();
        this->ComputePartTopology();
        if (usingOCC)
        {
            this->ComputePartGeometry();
            //status &= this->ComputeCurveSurfaceGeometryMaps();
        }
    }
    else
        status = false;
    return status;
}


int OCC_Computation::ReadFile(const char* cadFilePath)
{
    auto ClearWorkSession = [](const Handle(XSControl_WorkSession)& workSession) -> void
    {
        if (workSession.IsNull()) return;

        // clear the transient process
        const Handle(Transfer_TransientProcess)& mapReader = workSession->TransferReader()->TransientProcess();
        if (!mapReader.IsNull())
            mapReader->Clear();

        // clear the transfer reader
        const Handle(XSControl_TransferReader)& transferReader = workSession->TransferReader();
        if (!transferReader.IsNull())
            transferReader->Clear(-1);
    };

    // [&] captures *this and all local variables
    auto ReadOperation = [&](   XSControl_Reader* reader,
                                bool igesFormat ) -> int
    {
        Handle(XSControl_WorkSession) workSession = reader->WS();
        IFSelect_ReturnStatus status = reader->ReadFile(cadFilePath);
        if (IFSelect_RetDone != status) return -1;
        int numTranslations = reader->TransferRoots();
        if (!numTranslations) return -1;
        auto preResult = reader->OneShape();
        if (preResult.IsNull()) return -1;
        workSession->NewModel();
        Standard::Purge();
        *this->m_Shape = std::move(preResult);
        if (igesFormat)
            ClearWorkSession(workSession);
        return 0;
    };

    if (!cadFilePath)
        cadFilePath = this->m_FilePath.c_str();
    int status = 0;
    bool igesFormat = false;
    auto fileExtn = Utility::FileExtn(cadFilePath);
    if ("step" == fileExtn || "stp" == fileExtn)
    {
        STEPControl_Reader stpReader;
        status = ReadOperation(&stpReader, igesFormat);
    }
    else if ("iges" == fileExtn || "igs" == fileExtn)
    {
        IGESControl_Reader igsReader;
        status = ReadOperation(&igsReader, igesFormat = true);
    }
    else if ("brep" == fileExtn)
    {
        BRep_Builder builder;
        status &= BRepTools::Read(*this->m_Shape, cadFilePath, builder) ? 0 : -1;
    }
    else
    {
        return -1;
    }

    if (status) return -1;
#if 0
    IFSelect_PrintCount mode = IFSelect_ListByItem;
    reader.PrintCheckLoad(false, mode);
    Standard_Integer numRoots = reader.NbRootsForTransfer();
    Standard_Integer numTranslations = reader.TransferRoots();
    std::cout << numRoots << std::endl;
    std::cout << numTranslations << std::endl;
    BRepTools::Write(*this->m_Shape, "output.brep");
    //TopoDS_Shape shape = reader.Shape();
#endif

    if (igesFormat)
    {
        // topo-geom maps are needed for computation of sewing tolerance
        status &= this->GenerateTopologyGeometryMaps() ? 0 : -1;
        status &= this->ConvertIGESToSolid() ? 0 : -1;
        if (Utility::DEBUG_LEVEL >= 1)
            std::cout << "\n\nIGES ShapeType: " << \
            TopAbs::ShapeTypeToString(this->m_Shape->ShapeType()) << std::endl;
    }

    return status;
}


int OCC_Computation::WriteFile(int16_t format)
{
    int status = 0;
    if (this->m_FilePath.empty() || this->m_Shape->IsNull())
        return -1;

    auto filePath = this->m_FilePath.c_str();
    auto fileExtn = Utility::FileExtn(filePath);
    // prohibit file writing if extension is identical
    if (((0 == format) && ("step" == fileExtn || "stp" == fileExtn)) ||
        ((1 == format) && ("iges" == fileExtn || "igs" == fileExtn)) ||
        ((2 == format) && "brep" == fileExtn))
        return -1;
    auto outputPath = Utility::OutputPath(filePath, format);

    if (0 == format)
    {
        STEPControl_Writer stpWriter;
        IFSelect_ReturnStatus status = stpWriter.Transfer(*this->m_Shape, STEPControl_ManifoldSolidBrep);
        if (IFSelect_RetDone != status)
            status = stpWriter.Transfer(*this->m_Shape, STEPControl_BrepWithVoids);
        if (IFSelect_RetDone != status) return -1;
        status = stpWriter.Write(outputPath.c_str());
        if (IFSelect_RetDone != status) return -1;
    }
    else if (1 == format)
    {
        IGESControl_Writer igsWriter;
        bool status = igsWriter.AddShape(*this->m_Shape);
        if (!status) return -1;
        status = igsWriter.Write(outputPath.c_str());
        if (!status) return -1;
    }
    else if (2 == format)
    {
        status &= BRepTools::Write(*this->m_Shape, outputPath.c_str()) ? 0 : -1;
    }
    else
        return -1;

    return 0;
}


bool OCC_Computation::ConvertIGESToSolid()
{
    auto MinMaxAvgeTolerance = [](  const std::shared_ptr<TopTools_IndexedMapOfShape>& topoMap,
                                    TopAbs_ShapeEnum topoType, double& min, double& max, double& avge) -> void
    {
        for (int idx = 1; idx <= topoMap->Extent(); ++idx)
        {
            auto ithEntity = (*topoMap)(idx);
            auto ithTolerance = 0.0;
            if (TopAbs_VERTEX == topoType)
                ithTolerance = BRep_Tool::Tolerance(TopoDS::Vertex(ithEntity));
            else
                ithTolerance = BRep_Tool::Tolerance(TopoDS::Edge(ithEntity));
            if (ithTolerance < min)
                min = ithTolerance;
            if (ithTolerance > max)
                max = ithTolerance;
            //std::cout << ithTolerance << std::endl;
            avge += ithTolerance;
        }
        avge /= topoMap->Extent();
    };

    // return if the shape is already a solid or if topology maps are not generated
    if (TopAbs_SOLID == this->m_Shape->ShapeType() || 
        this->m_VertexMap->IsEmpty() || this->m_EdgeMap->IsEmpty())
       return false;

#if 1
    double minVertexTolerance = DBL_MAX, minEdgeTolerance = DBL_MAX, avgeVertexTolerance = 0.0;
    double maxVertexTolerance = DBL_MIN, maxEdgeTolerance = DBL_MIN, avgeEdgeTolerance = 0.0;
    MinMaxAvgeTolerance(this->m_VertexMap, TopAbs_VERTEX, minVertexTolerance, maxVertexTolerance, avgeVertexTolerance);
    MinMaxAvgeTolerance(this->m_EdgeMap, TopAbs_EDGE, minEdgeTolerance, maxEdgeTolerance, avgeEdgeTolerance);
    auto maxFaceTolerance = BRep_Tool::MaxTolerance(*this->m_Shape, TopAbs_FACE);

    auto minTol = std::min(avgeVertexTolerance, avgeEdgeTolerance);
    auto maxTol = std::max(1e1 * maxVertexTolerance, 1e1 * maxEdgeTolerance);
    auto boundingBoxDiagonal = this->ComputeBoundingBox(0.0f, true).Modulus();
    // criteria can be modified
    // the tolerance is set slightly on the higher end to ensure that the resultant sewed shape
    // consists of a single shell
    this->m_Tolerance = std::min(maxTol, 1e-3 * boundingBoxDiagonal);
#else
    ShapeAnalysis_ShapeTolerance shapeTolerance[2];
    shapeTolerance[0].InitTolerance();
    shapeTolerance[1].InitTolerance();
    // 3 => min, avge, max
    double vertexTol[3]{ 0.0 }, edgeTol[3]{ 0.0 };
    shapeTolerance[0].AddTolerance(*this->m_Shape, TopAbs_VERTEX);
    shapeTolerance[1].AddTolerance(*this->m_Shape, TopAbs_EDGE);
    // -1 => min, = => avge, 1 => max
    for (int8_t i = -1, j = 0; i <= 1; ++i, ++j)
    {
        vertexTol[j] = shapeTolerance[0].GlobalTolerance(i);
        edgeTol[j] = shapeTolerance[1].GlobalTolerance(i);
    }
    auto minTol = std::min(vertexTol[0], edgeTol[0]);
    auto maxTol = std::max(vertexTol[2], edgeTol[2]);
    auto avgeTol = 0.5 * (vertexTol[1] + edgeTol[1]);
    this->m_Tolerance = std::max(vertexTol[1], edgeTol[1]);
#endif

    BRepBuilderAPI_Sewing sewFaces(this->m_Tolerance);
    sewFaces.SetMinTolerance(minTol);
    sewFaces.SetMaxTolerance(maxTol);
    // combine multiple faces for BRepBuilderAPI_Sewing
    for (TopExp_Explorer explorer(*this->m_Shape, TopAbs_FACE); explorer.More(); explorer.Next())
    {
        if (!explorer.Current().IsNull())
            sewFaces.Add(explorer.Current());
    }
    sewFaces.Perform();
    const TopoDS_Shape& sewedShape = sewFaces.SewedShape();
    if (sewedShape.IsNull()) return false;

    if (TopAbs_COMPOUND == sewedShape.ShapeType())
    {
        // sewedShape consists of multiple shells - this can happen if stitching tolerance is less
        TopoDS_Shell shell;
        BRep_Builder shellBuilder;
        shellBuilder.MakeShell(shell);
        for (TopExp_Explorer explorer(*this->m_Shape, TopAbs_FACE); explorer.More(); explorer.Next())
        {
            if (!explorer.Current().IsNull())
                shellBuilder.Add(shell, explorer.Current());
        }

        BRepBuilderAPI_Sewing sewFacesFromShells(1e-4);
        sewFacesFromShells.Add(shell);
        sewFacesFromShells.Perform();
        TopoDS_Shape sewedShapeFromShells = sewFacesFromShells.SewedShape();

        std::vector<TopoDS_Shell> shells;
        for (TopExp_Explorer explorer(sewedShapeFromShells, TopAbs_SHELL); explorer.More(); explorer.Next())
        {
            if (!explorer.Current().IsNull())
                shells.emplace_back(TopoDS::Shell(explorer.Current()));
        }

        TopoDS_Solid solid;
        BRep_Builder solidBuilder;
        solidBuilder.MakeSolid(solid);
        for (const auto& itr : shells)
            solidBuilder.Add(solid, itr);

        BRepCheck_Analyzer sanityCheck(solid, true, true);
        if (!sanityCheck.IsValid())
            return false;
        this->m_Shape = std::make_shared<TopoDS_Shape>(solid);
    }
    else if (TopAbs_SHELL == sewedShape.ShapeType())
    {
        // sewedShape consists of a single shell
        BRepBuilderAPI_MakeSolid solid;
        solid.Add(TopoDS::Shell(sewedShape));
        this->m_Shape = std::make_shared<TopoDS_Shape>(solid.Shape());
    }
    else
        return false;

    return true;
}


std::pair<int16_t, int16_t> OCC_Computation::CountProceduralGeometryEntities()
{
    int16_t count = 0;
    std::pair<int16_t, int16_t> proceduralEntities;
    auto citr = this->m_CurveGeometryMap->find(GeomAbs_CurveType::GeomAbs_OffsetCurve);
    if (this->m_CurveGeometryMap->end() != citr)
        proceduralEntities.first = citr->second;
    std::for_each(  this->m_SurfaceGeometryMap->begin(), this->m_SurfaceGeometryMap->end(),
                    [&count](auto surfaceType) -> int16_t
                    {
                        if (GeomAbs_SurfaceType::GeomAbs_SurfaceOfRevolution == surfaceType.first)
                            count += surfaceType.second;
                        else if (GeomAbs_SurfaceType::GeomAbs_SurfaceOfExtrusion == surfaceType.first)
                            count += surfaceType.second;
                        else if (GeomAbs_SurfaceType::GeomAbs_OffsetSurface == surfaceType.first)
                            count += surfaceType.second;
                        return count;
                    });
    proceduralEntities.second = count;
    return proceduralEntities;
}


std::pair<int16_t, int16_t> OCC_Computation::ConvertProceduralGeometryToBSplineGeometry()
{
    if (!this->m_TopoGeomMapsGenerated) return { -1, -1 };
    std::pair<int16_t, int16_t> convertedEntities;

    // procedural curves : offset
    for (int16_t idx = 1; idx <= this->m_EdgeMap->Extent(); ++idx)
    {
        const TopoDS_Edge& edge = TopoDS::Edge((*this->m_EdgeMap)(idx));
        BRepAdaptor_Curve brepCurve(edge);

        if (GeomAbs_OffsetCurve == brepCurve.GetType())
        {
            auto tMin = 0.0, tMax = 0.0;
            Handle(Geom_Curve) geomCurve = BRep_Tool::Curve(edge, tMin, tMax);
            if (geomCurve.IsNull())
                continue;
            Handle(Geom_TrimmedCurve) trimmedCurve = new Geom_TrimmedCurve(geomCurve, tMin, tMax);
            Handle(Geom_BSplineCurve) bsplineCurve = GeomConvert::CurveToBSplineCurve(trimmedCurve);
            // TBD: needs additional changes for rebuilding geometry/topology
            BRep_Builder edgeBuilder;
            edgeBuilder.UpdateEdge(edge, bsplineCurve, edge.Location(), BRep_Tool::Tolerance(edge));
            ++convertedEntities.first;
        }
    }

    // procedural curves : offset, linear-extrusion, revolution
    for (int16_t idx = 1; idx <= this->m_FaceMap->Extent(); ++idx)
    {
        const TopoDS_Face& face = TopoDS::Face((*this->m_FaceMap)(idx));
        BRepAdaptor_Surface brepSurface(face);
        auto type = brepSurface.GetType();

        if (GeomAbs_SurfaceOfRevolution == type ||
            GeomAbs_SurfaceOfExtrusion == type ||
            GeomAbs_OffsetSurface == type)
        {
            auto uMin = brepSurface.FirstUParameter();
            auto uMax = brepSurface.LastUParameter();
            auto vMin = brepSurface.FirstVParameter();
            auto vMax = brepSurface.LastVParameter();
            // the below code is computationally expensive, hence commented out
            //double uMin = 0.0, uMax = 0.0, vMin = 0.0, vMax = 0.0;
            //BRepTools::UVBounds(face, uMin, uMax, vMin, vMax);
            Handle(Geom_Surface) geomSurface = BRep_Tool::Surface(face);
            if (geomSurface.IsNull())
                continue;
            Handle(Geom_RectangularTrimmedSurface) trimmedSurface =
                new Geom_RectangularTrimmedSurface(geomSurface, uMin, uMax, vMin, vMax);
            Handle(Geom_BSplineSurface) bsplineSurface = GeomConvert::SurfaceToBSplineSurface(trimmedSurface);
            // TBD: needs additional changes for rebuilding geometry/topology
            BRep_Builder faceBuilder;
            faceBuilder.UpdateFace(face, bsplineSurface, face.Location(), BRep_Tool::Tolerance(face));
            ++convertedEntities.second;
        }
    }

    BRepCheck_Analyzer sanityCheck(*this->m_Shape, true, true);
    if (!sanityCheck.IsValid()) return { -1, -1 };
    return convertedEntities;
}


std::array<int16_t, 3> OCC_Computation::FixRedundantTopology()
{
    if (!this->m_TopoGeomMapsGenerated) return {-1, -1, -1};
    std::array<int16_t, 3> reducedEntities{ 0 };

    // the tolerances are intentionally loose because canonical conversion of IGES CAD files
    // can introduce some tolerance errors
    constexpr double linTol = 1e-1, angTol = 1e-2;
    auto maxCumulativeEdgeLength = 0.0f;

    // edge group data : radius, axis, shared faces, cumulative edge length
    // edges which are to be combined are grouped together
    // cumulative edge length is the aggregate length of individual tiny edges
    // that can be combined together and it is a parameter which is to be set
    // programmatically to merge edges. ideally each edge group should have a
    // distinct value and should be merged accordingly.
    std::vector<std::tuple<float, gp_Ax1, int16_t, int16_t, float>> groupedEdges;
    for (const auto& itr : *this->m_EdgeTypeMap)
    {
        bool newEdgeGroup = true;
        const TopoDS_Edge& edge = TopoDS::Edge((*this->m_EdgeMap)(itr.first));
        BRepAdaptor_Curve curve(edge);

        auto sharedFaces = this->m_EdgeFaceMap->FindFromIndex(itr.first);
        auto thisFace1 = this->m_FaceMap->FindIndex(sharedFaces.First());
        auto thisFace2 = this->m_FaceMap->FindIndex(sharedFaces.Last());

        GProp_GProps edgeProps;
        BRepGProp::LinearProperties(edge, edgeProps);
        auto thisEdgeLength = edgeProps.Mass();

        gp_Ax1 thisAxis;
        gp_Dir thisDirection;
        gp_Pnt thisLocation;
        float thisRadius = 0.0f;

        // line & circle are supported for now
        if (GeomAbs_CurveType::GeomAbs_Line == itr.second)
        {
            gp_Lin line = curve.Line();
            thisAxis.SetDirection(line.Direction());
            thisAxis.SetLocation(line.Location());
            thisDirection = line.Direction();
            thisLocation = line.Location();
        }
        else if (GeomAbs_CurveType::GeomAbs_Circle == itr.second)
        {
            gp_Circ circle = curve.Circle();
            thisAxis = circle.Axis();
            thisRadius = static_cast<float>(circle.Radius());
        }

        for (auto& [rad, ax, fc1, fc2, len] : groupedEdges)
        {
            // edges to be combined should necessarily have identical shared faces
            if (fc1 == thisFace1 && fc2 == thisFace2)
            {
                if (GeomAbs_CurveType::GeomAbs_Line == itr.second)
                {
                    // linear edges can be combined only if each line's axis (location & direction)
                    // is identical
                    // (ax.Location().Distance(thisLocation) < distTol)
                    newEdgeGroup = !(ax.Direction().IsParallel(thisDirection, angTol) &&
                        gp_Vec(ax.Location(), thisLocation).IsParallel(gp_Vec(thisAxis.Direction()), angTol));
                }
                else if (GeomAbs_CurveType::GeomAbs_Circle == itr.second)
                {
                    // circular edges can be combined only if each circles's axis (center & normal) and radius
                    // is identical
                    newEdgeGroup = !((fabs(rad - thisRadius) < linTol) &&
                        ax.IsCoaxial(thisAxis, angTol, linTol));
                }
                if (!newEdgeGroup)
                {
                    len += static_cast<float>(thisEdgeLength);
                    if (len > maxCumulativeEdgeLength)
                        maxCumulativeEdgeLength = len;
                    break;
                }
            }
        }
        if (newEdgeGroup)
        {
            // std::make_tuple()
            groupedEdges.emplace_back(static_cast<float>(thisRadius),
                thisAxis, thisFace1, thisFace2, static_cast<float>(thisEdgeLength));
        }
    }

    if (maxCumulativeEdgeLength > Precision::Confusion())
    {
        ShapeFix_Wireframe mergeEdges(*this->m_Shape);
        // this value determines edge geometry (ShapeFix_Wireframe::MergeSmallEdges) and consequently those
        // edges which are eligible for merging
        // angular tolerance should be greater than:
        // at each shared vertex angle between tangents of successive edges (along the curve)
        mergeEdges.SetLimitAngle(angTol);                         // 0.01 rad. = 0.57 deg.
        // this value determines small edges (ShapeAnalysis_Wire::CheckSmall) and consequently those
        // edges which are eligible for merging
        // distance tolerance should be greater than:
        // straight line distance between 2 vertices of short edge (curve geometry ignored)
        // straight line distance between each vertex to edge mid-point (along the curve)
        mergeEdges.SetPrecision(maxCumulativeEdgeLength);
        mergeEdges.SetContext(new ShapeBuild_ReShape);
        if (mergeEdges.FixSmallEdges())
        {
            TopoDS_Shape modifiedShape = mergeEdges.Shape();
            BRepCheck_Analyzer sanityCheck(modifiedShape, true, true);
            if (!sanityCheck.IsValid()) return { -1, -1, -1 };

            this->m_Shape = std::make_shared<TopoDS_Shape>(modifiedShape);
            ShapeAnalysis_ShapeContents shapeCount;
            shapeCount.Perform(*this->m_Shape);
            reducedEntities[0] = this->m_VertexMap->Extent() - shapeCount.NbSharedVertices();
            reducedEntities[1] = this->m_EdgeMap->Extent() - shapeCount.NbSharedEdges();
            reducedEntities[2] = 0;
        }
    }

    // TBD: merging of faces

    return reducedEntities;
}


gp_XYZ OCC_Computation::ComputeBoundingBox( float offset,
                                            bool igesMode   )
{
    auto ComputeAxesDim = [](   const Bnd_Box& aabb,
                                float (&dimensions)[3]  ) -> void
    {
        dimensions[0] = static_cast<float>(fabs(aabb.CornerMin().X() - aabb.CornerMax().X()));
        dimensions[1] = static_cast<float>(fabs(aabb.CornerMin().Y() - aabb.CornerMax().Y()));
        dimensions[2] = static_cast<float>(fabs(aabb.CornerMin().Z() - aabb.CornerMax().Z()));
    };

#if 0
    GProp_GProps propertiesSystemFace;
    BRepGProp::VolumeProperties(shape, propertiesSystemFace);
    double shapeVolume = propertiesSystemFace.Mass();
    gp_Pnt centerMass = propertiesSystemFace.CentreOfMass();
    std::cout << shapeVolume << std::endl;
    std::cout << "(" << centerMass.X() << ", " << centerMass.Y() << ", " << centerMass.Z() << ")" << std::endl;
#endif

    float axesDims[3]{ 0.0f }, pcaDims[3]{ 0.0f }, offAxesDims[3]{ 0.0f };
    if (!igesMode || this->m_AxisAlignedBoundingBox->IsVoid())
    {
        // compute AABB only once when for IGES parts
        BRepBndLib::Add(*this->m_Shape, *this->m_AxisAlignedBoundingBox);
        ComputeAxesDim(*this->m_AxisAlignedBoundingBox, axesDims);
    }

    gp_XYZ offsettedBoundingBoxDim(axesDims[0], axesDims[1], axesDims[2]);
    if (!igesMode)
    {
        // query AABB dimensions before enlarging the box
        ComputeAxesDim(*this->m_AxisAlignedBoundingBox, axesDims);

        // PCA is used for finding PMI (principal moment of inertia) axes
        BRepBndLib::AddOBB(*this->m_Shape, *this->m_OrientedBoundingBox, false, true, true);
        // OBB generates box center and half dimensions along each axes
        pcaDims[0] = 2.0f * static_cast<float>(this->m_OrientedBoundingBox->XHSize());
        pcaDims[1] = 2.0f * static_cast<float>(this->m_OrientedBoundingBox->YHSize());
        pcaDims[2] = 2.0f * static_cast<float>(this->m_OrientedBoundingBox->ZHSize());

        // offsetted bounding box
        Bnd_Box offsettedBoundingBox;
        offsettedBoundingBox.Add(*this->m_AxisAlignedBoundingBox);
        offsettedBoundingBox.Enlarge(offset);
        ComputeAxesDim(offsettedBoundingBox, offAxesDims);
        offsettedBoundingBoxDim.SetX(offAxesDims[0]);
        offsettedBoundingBoxDim.SetY(offAxesDims[1]);
        offsettedBoundingBoxDim.SetZ(offAxesDims[2]);

        if (Utility::DEBUG_LEVEL >= 1)
        {
            std::cout << "Axis-Aligned-Bounding-Box: (" << axesDims[0] << ", " << axesDims[1] << ", " << axesDims[2] << ")" << std::endl;
            std::cout << "Oriented-Bounding-Box: (" << pcaDims[0] << ", " << pcaDims[1] << ", " << pcaDims[2] << ")" << std::endl;
            //std::cout << "OBB == AABB: " << std::boolalpha << this->m_OrientedBoundingBox->IsAABox() << std::endl;
            std::cout << "Offsetted-Bounding-Box: (" << \
                offsettedBoundingBoxDim.X() << ", " << offsettedBoundingBoxDim.Y() << ", " << offsettedBoundingBoxDim.Z() << ")" << std::endl;
        }
    }

    return offsettedBoundingBoxDim;
}


void OCC_Computation::ComputePartGeometry(bool usingOCC) const
{
    if (this->m_VertexMap->IsEmpty()) return;

    if (usingOCC)
        std::cout << "\nVertex-ID-Point Map:" << std::endl;
    for (int16_t idx = 1; idx <= this->m_VertexMap->Extent(); ++idx)
    {
        const TopoDS_Vertex& vertex = TopoDS::Vertex((*this->m_VertexMap)(idx));
        gp_Pnt point = BRep_Tool::Pnt(vertex);
        if (usingOCC)
        {
            std::cout << std::fixed << std::setprecision(2) << idx << \
                ": Coordinates: (" << point.X() << ", " << point.Y() << ", " << point.Z() << ")" << std::endl;
        }
    }

    int16_t numDegenerateEdges = 0;
    if (usingOCC)
        std::cout << "\nEdge-ID-Geometry Map:" << std::endl;
    for (int16_t idx = 1; idx <= this->m_EdgeMap->Extent(); ++idx)
    {
        const TopoDS_Edge& edge = TopoDS::Edge((*this->m_EdgeMap)(idx));
        bool isDegenerate = false;
        if (BRep_Tool::Degenerated(edge))
        {
            isDegenerate = true;
            ++numDegenerateEdges;
        }

        BRepAdaptor_Curve curve(edge);
        ShapeAnalysis_Edge edgeAnalysis;
        bool isClosed = edgeAnalysis.IsClosed3d(edge);
        Handle(Geom_Curve) curve3d;
        double tStart = 0.0, tLast = 0.0, tPeriod = 0.0, length = 0.0;
        bool isPeriodic = false, isSeam = false;
        if (edgeAnalysis.Curve3d(edge, curve3d, tStart, tLast))
            isPeriodic = ShapeAnalysis_Curve::IsPeriodic(curve3d);
        if (isPeriodic)
            tPeriod = curve3d->Period();
        if (!isDegenerate)
        {
            // seams can also be computed using ShapeExtend_WireData class
            const TopoDS_Face& face = TopoDS::Face(this->m_EdgeFaceMap->FindFromKey(edge).First());
            isSeam = edgeAnalysis.IsSeam(edge, face);
            GProp_GProps edgeProps;
            BRepGProp::LinearProperties(edge, edgeProps);
            length = edgeProps.Mass();
        }
        auto crvType = OCC_Computation::CurveGeometry(curve.GetType());
        if (usingOCC)
            std::cout << idx << ": Geometry: " << crvType;
        if (isDegenerate)
        {
            gp_Pnt point = curve.Value(curve.FirstParameter());
            if (usingOCC)
            {
                std::cout << ", Degenerate: true, ";
                std::cout << std::fixed << std::setprecision(2) << \
                    "Coordinates: (" << point.X() << ", " << point.Y() << ", " << point.Z() << ")";
            }
        }
        else
        {
            if (usingOCC)
            {
                std::cout << ", Length = " << length;
                std::cout << std::boolalpha << ", Closed: " << isClosed;
                if (isSeam)
                    std::cout << ", Seam: true";
                if (isPeriodic)
                    std::cout << ", Periodic: true, Period = " << tPeriod;
            }
        }
        if (usingOCC) std::cout << std::endl;
    }

    if (usingOCC)
        std::cout << "\nFace-ID-Geometry Map:" << std::endl;
    for (int16_t idx = 1; idx <= this->m_FaceMap->Extent(); ++idx)
    {
        TopTools_IndexedMapOfShape verticesOnFaceMap;
        const TopoDS_Face& face = TopoDS::Face((*this->m_FaceMap)(idx));
        BRepAdaptor_Surface surface(face);
        int8_t numVertices = 0;
        for (TopExp_Explorer exp(face, TopAbs_VERTEX); exp.More(); exp.Next())
        {
            const TopoDS_Vertex& vertex = TopoDS::Vertex(exp.Current());
            verticesOnFaceMap.Add(vertex);
        }
        double uPeriod = 0.0, vPeriod = 0.0, area = 0.0;
        bool isUClosed = surface.IsUClosed();
        bool isVClosed = surface.IsVClosed();
        bool isUPeriodic = surface.IsUPeriodic();
        bool isVPeriodic = surface.IsVPeriodic();
        if (isUPeriodic)
            uPeriod = surface.UPeriod();
        if (isVPeriodic)
            vPeriod = surface.VPeriod();
        auto srfType = OCC_Computation::SurfaceGeometry(surface.GetType());
        GProp_GProps faceProps;
        BRepGProp::SurfaceProperties(face, faceProps);
        area = faceProps.Mass();
        if (usingOCC)
        {
            std::cout << idx << ": Geometry: " << srfType << ", #Vertices = " << verticesOnFaceMap.Extent() << \
                ", Area = " << area << std::boolalpha << ", U-Closed: " << isUClosed << ", V-Closed: " << isVClosed;
            if (isUPeriodic)
                std::cout << ", U-Periodic: true, U-Period = " << uPeriod;
            if (isVPeriodic)
                std::cout << ", V-Periodic: true, V-Period = " << vPeriod;
            std::cout << std::endl;
        }
    }
}


bool OCC_Computation::ComputeCurveSurfaceGeometryMaps()
{
    auto UpdateMaps = [&](  GeomAbs_CurveType curveType,
                            GeomAbs_SurfaceType surfaceType ) -> void
    {
        if (GeomAbs_SurfaceType::GeomAbs_OtherSurface == surfaceType)
        {
            auto itr = this->m_CurveGeometryMap->find(curveType);
            if (itr == this->m_CurveGeometryMap->end())
                this->m_CurveGeometryMap->emplace(std::make_pair(curveType, 1));
            else
                itr->second += 1;
        }
        else if (GeomAbs_CurveType::GeomAbs_OtherCurve == curveType)
        {
            auto itr = this->m_SurfaceGeometryMap->find(surfaceType);
            if (itr == this->m_SurfaceGeometryMap->end())
                this->m_SurfaceGeometryMap->emplace(std::make_pair(surfaceType, 1));
            else
                itr->second += 1;
        }
    };

    if (this->m_VertexMap->IsEmpty()) return false;
    if (!this->m_CurveGeometryMap->empty())
        this->m_CurveGeometryMap->clear();
    if (!this->m_SurfaceGeometryMap->empty())
        this->m_SurfaceGeometryMap->clear();

    for (int16_t idx = 1; idx <= this->m_EdgeMap->Extent(); ++idx)
    {
        int16_t count = 0;
        const TopoDS_Edge& edge = TopoDS::Edge((*this->m_EdgeMap)(idx));
        BRepAdaptor_Curve curve(edge);
        if (GeomAbs_CurveType::GeomAbs_Line == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_Line, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_Circle == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_Circle, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_Ellipse == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_Ellipse, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_Hyperbola == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_Hyperbola, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_Parabola == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_Parabola, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_BezierCurve == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_BezierCurve, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_BSplineCurve == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_BSplineCurve, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
        else if (GeomAbs_CurveType::GeomAbs_OffsetCurve == curve.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OffsetCurve, GeomAbs_SurfaceType::GeomAbs_OtherSurface);
    }

    for (int16_t idx = 1; idx <= this->m_FaceMap->Extent(); ++idx)
    {
        const TopoDS_Face& face = TopoDS::Face((*this->m_FaceMap)(idx));
        BRepAdaptor_Surface surface(face);
        if (GeomAbs_SurfaceType::GeomAbs_Plane == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_Plane);
        else if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_Cylinder);
        else if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_Cone);
        else if (GeomAbs_SurfaceType::GeomAbs_Sphere == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_Sphere);
        else if (GeomAbs_SurfaceType::GeomAbs_Torus == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_Torus);
        else if (GeomAbs_SurfaceType::GeomAbs_BezierSurface == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_BezierSurface);
        else if (GeomAbs_SurfaceType::GeomAbs_BSplineSurface == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_BSplineSurface);
        else if (GeomAbs_SurfaceType::GeomAbs_SurfaceOfRevolution == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_SurfaceOfRevolution);
        else if (GeomAbs_SurfaceType::GeomAbs_SurfaceOfExtrusion == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_SurfaceOfExtrusion);
        else if (GeomAbs_SurfaceType::GeomAbs_OffsetSurface == surface.GetType())
            UpdateMaps(GeomAbs_CurveType::GeomAbs_OtherCurve, GeomAbs_SurfaceType::GeomAbs_OffsetSurface);
    }

    return true;
}


std::pair<int16_t, int16_t> OCC_Computation::ComputePartTopology()
{
    if (this->m_EdgeMap->IsEmpty()) return { -1, -1 };

    int16_t numDegenerateEdges = 0, numSeams = 0;
    ShapeAnalysis_Edge edgeAnalysis;
    for (int16_t idx = 1; idx <= this->m_EdgeMap->Extent(); ++idx)
    {
        const TopoDS_Edge& edge = TopoDS::Edge((*this->m_EdgeMap)(idx));
        if (BRep_Tool::Degenerated(edge))
        {
            // may or may not be a pole
            ++numDegenerateEdges;
        }
        else
        {
            const TopoDS_Face& face = TopoDS::Face(this->m_EdgeFaceMap->FindFromKey(edge).First());
            if (edgeAnalysis.IsSeam(edge, face))
                ++numSeams;
        }
    }

    TopTools_ShapeSet topology;
    topology.Add(*this->m_Shape);
    if (Utility::DEBUG_LEVEL >= 1)
    {
        std::cout << "\nCAD Topology:" << std::endl;
        topology.DumpExtent(std::cout);
    }
    // Euler-Poincare formula: V - E + F - (L - F) - 2*(S - G) = 0
    // multiloops = M = [L - F]
    // genus = G = [S - (V - E + F - M) / 2]
    ShapeAnalysis_ShapeContents shapeCount;
    shapeCount.Perform(*this->m_Shape);
    //shapeCount.NbFaceWithSevWires();      // number of faces containing multiple loops
    // faces containing multiple loops may or may not contain holes
    decltype(numDegenerateEdges) numMultiLoops = shapeCount.NbSharedWires() - shapeCount.NbSharedFaces();
    this->m_Genus = shapeCount.NbSharedShells() - \
        (shapeCount.NbSharedVertices() - (shapeCount.NbSharedEdges() - numDegenerateEdges) + shapeCount.NbSharedFaces() - numMultiLoops) / 2;
    // factor 2 is used because single genus comprises of 2 holes for the CAD parts being used
    this->m_Holes = this->m_Genus > 0 ? (numMultiLoops - (2 * this->m_Genus)) : numMultiLoops;
    if (Utility::DEBUG_LEVEL >= 1)
    {
        std::cout << "DEGENERECIES: " << numDegenerateEdges << ", MULTI-LOOPS = " << numMultiLoops << \
            ", HOLES <= " << this->m_Holes << ", GENUS = " << this->m_Genus << "\n\n";
    }
    return { numDegenerateEdges, numMultiLoops };
}


const char* OCC_Computation::CurveGeometry(GeomAbs_CurveType curveType)
{
    switch (curveType)
    {
    case 0:
        return "Line";
    case 1:
        return  "Circle";
    case 2:
        return  "Ellipse";
    case 3:
        return  "Hyperbola";
    case 4:
        return  "Parabola";
    case 5:
        return  "Bezier";
    case 6:
        return  "BSpline";
    case 7:
        return  "Offset";
    default:
        return  "Other";
    }
}


const char* OCC_Computation::SurfaceGeometry(GeomAbs_SurfaceType surfaceType)
{
    switch (surfaceType)
    {
    case 0:
        return "Plane";
    case 1:
        return "Cylinder";
    case 2:
        return "Cone";
    case 3:
        return "Sphere";
    case 4:
        return "Torus";
    case 5:
        return "Bezier";
    case 6:
        return "BSpline";
    case 7:
        return "Revolution";
    case 8:
        return "Extrusion";
    case 9:
        return "Offset";
    default:
        return "Other";
    }
}


gp_Pnt OCC_Computation::ComputeCenter(const BRepAdaptor_Surface& surface)
{
    auto uFirst = surface.FirstUParameter();
    auto uLast = surface.LastUParameter();
    auto vFirst = surface.FirstVParameter();
    auto vLast = surface.LastVParameter();

    gp_Ax1 cylinderAxis;
    if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
        cylinderAxis = surface.Cylinder().Axis();
    else if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
        cylinderAxis = surface.Cone().Axis();

    // center is on the axis, not on the curved surface of cylinder/cone
    gp_Pnt center;
    gp_Pnt midHeightPoint = surface.Value(0.5 * (uFirst + uLast), 0.5 * (vFirst + vLast));
    Handle(Geom_Line) axisLine = new Geom_Line(cylinderAxis);
    GeomAPI_ProjectPointOnCurve pointProjector(midHeightPoint, axisLine);
    center.SetXYZ(pointProjector.NearestPoint().XYZ());

    return center;
}


void OCC_Computation::ComputeDegenerateQuadricDimensions(   const BRepAdaptor_Surface& surface,
                                                            ConeCylinderDimensions* dimensions  )
{
    if (!dimensions) return;

    auto uFirst = surface.FirstUParameter();
    auto uLast = surface.LastUParameter();
    auto vFirst = surface.FirstVParameter();
    auto vLast = surface.LastVParameter();
    bool semiArc = true;
    double arcAngle = 0.0;

    // surface vector on the cylinder surface is formed such that it is parallel to the axis
    // base vector is formed such that it is along the base radius
    // top vector is necessary only when the surface is a frustum of a cone
    for (int16_t i = 0; i < 2; ++i)
    {
        gp_Vec surfVec, baseVec, topVec;
        if (surface.IsUPeriodic())
        {
            // u is along circular direction, v is along height
            if (fabs((uLast - uFirst) - 2 * M_PI) < Precision::Confusion())
            {
                // full cylinder of 2*pi radians in u direction
                // interval is halved, else it will give the same point in R3 space
                if (uFirst > Precision::Confusion())
                    uFirst *= 0.5;
                else if (uLast > Precision::Confusion())
                    uLast *= 0.5;
            }
            else if (fabs((uLast - uFirst) - M_PI) > Precision::Confusion())
            {
                // neither half cylinder, nor full cylinder in u direction
                semiArc = false;
                arcAngle = static_cast<float>(fabs(uLast - uFirst));
                // 0 < arcAngle < 180
                if (arcAngle > M_PI) arcAngle = 2 * M_PI - arcAngle;
            }

            if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
            {
                if (!i)
                {
                    surfVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uFirst, vLast)));
                    baseVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uLast, vFirst)));
                }
                else
                {
                    surfVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                    baseVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uLast, vLast)));
                }
            }
            else if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
            {
                if (surface.Cone().SemiAngle() > 0)
                {
                    if (vFirst < 0.0 && vLast < 1e3 * Precision::Confusion())
                    {
                        // vLast = 0 (circular-base), vFirst < 0 (apex)
                        if (!i)
                            surfVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uFirst, vFirst)));
                        else
                            surfVec.Add(gp_Vec(surface.Value(uLast, vLast), surface.Value(uLast, vFirst)));
                        baseVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uLast, vLast)));
                        topVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uLast, vFirst)));
                    }
                    else if (vFirst < 1e3 * Precision::Confusion() && vLast > 0.0)
                    {
                        // vFirst = 0 (circular-base), vLast > 0 (apex may or may not be present)
                        if (!i)
                            surfVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uFirst, vLast)));
                        else
                            surfVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                        baseVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uLast, vFirst)));
                        topVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uLast, vLast)));
                    }
                }
                else
                {
                    // vLast = 0 (circular-base), vFirst > 0 (apex) - to be implemented
                    // vFirst = 0 (circular-base), vLast > 0 (apex may or may not be present) - to be implemented
                    ;
                }
            }
        }
        else if (surface.IsVPeriodic())
        {
            // u is along height, v is along circular direction
            if (fabs((vLast - vFirst) - 2 * M_PI) < Precision::Confusion())
            {
                // full cylinder of 2*pi radians in v direction
                // interval is halved, else it will give the same point in R3 space
                if (vFirst > Precision::Confusion())
                    vFirst *= 0.5;
                else if (vLast > Precision::Confusion())
                    vLast *= 0.5;
            }
            else if (fabs((vLast - vFirst) - M_PI) > Precision::Confusion())
            {
                // neither half cylinder, nor full cylinder in u direction
                semiArc = false;
                dimensions->angle = static_cast<float>(fabs(vLast - vFirst));
                // 0 < angle < 180
                if (dimensions->angle > M_PI)
                    dimensions->angle = 2 * M_PI - dimensions->angle;
            }

            if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
            {
                if (!i)
                {
                    surfVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uLast, vFirst)));
                    baseVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uFirst, vLast)));
                }
                else
                {
                    surfVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uLast, vLast)));
                    baseVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                }
            }
            else if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
            {
                if (surface.Cone().SemiAngle() > 0)
                {
                    if (uFirst < 0.0 && uLast < 1e3 * Precision::Confusion())
                    {
                        // uLast = 0 (circular-base), uFirst < 0 (apex)
                        if (!i)
                            surfVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uFirst, vFirst)));
                        else
                            surfVec.Add(gp_Vec(surface.Value(uLast, vLast), surface.Value(uFirst, vLast)));
                        baseVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                        topVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uFirst, vLast)));
                    }
                    else if (uFirst < 1e3 * Precision::Confusion() && uLast > 0.0)
                    {
                        // uFirst = 0 (circular-base), uLast > 0 (apex may or may not be present)
                        if (!i)
                            surfVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uFirst, vLast)));
                        else
                            surfVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                        baseVec.Add(gp_Vec(surface.Value(uFirst, vFirst), surface.Value(uFirst, vLast)));
                        topVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                    }
                }
                else
                {
                    // uLast = 0 (circular-base), uFirst > 0 (apex) - to be implemented
                    // uFirst = 0 (circular-base), uLast != 0 (apex may or may not be present) - to be implemented
                    ;
                }
            }
        }

        // radius is computed from first surface uv parameterization instead of querying it directly
        // from the GeomSurface because the parameters don't always get updated after canonical
        // conversion (when IGES parts are provided)
        if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
        {
            gp_Cylinder cylSurface = surface.Cylinder();
            //dimensions->center = OCC_Computation::ComputeCenter(surface);
            // baseVec is checked for orthogonality against cylSurface.Axis() instead of
            // cylSurface.XAxis() / cylSurface.YAxis() because baseVec orientation can be
            // different and angular tolerance of 1e-10 is used to offset floating point errors
            // baseVec.IsNormal(cylSurface.XAxis().Direction(), 1e2 * Precision::Angular())
            // baseVec.IsNormal(cylSurface.YAxis().Direction(), 1e2 * Precision::Angular())
            if (surfVec.IsParallel(cylSurface.Axis().Direction(), 1e2 * Precision::Angular()) &&
                baseVec.IsNormal(cylSurface.Axis().Direction(), 1e2 * Precision::Angular()))
            {
                if (!semiArc)
                {
                    // cosine rule: computing radius from base vector
                    // r = vec / (sqrt(2 * (1 - cos(theta))))
                    auto num = baseVec.Magnitude();
                    auto den = std::sqrt(2 * (1 - std::cos(arcAngle)));
                    dimensions->radius = static_cast<float>(num / den);
                }
                else
                {
                    // normally CAD parts contain semi arc cylinders (base-angle = 180) and hence
                    // baseVec is along the diameter
                    dimensions->radius = static_cast<float>(0.5 * baseVec.Magnitude());
                }
                dimensions->height = static_cast<float>(surfVec.Magnitude());
                break;
            }
        }
        else if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
        {
            gp_Cone conSurface = surface.Cone();
            // baseVec is checked for orthogonality against conSurface.Axis() instead of
            // conSurface.XAxis() / conSurface.YAxis() because baseVec orientation can be
            // different and angular tolerance of 1e-10 is used to offset floating point errors
            // baseVec.IsNormal(conSurface.XAxis().Direction(), 1e2 * Precision::Angular())
            // baseVec.IsNormal(conSurface.YAxis().Direction(), 1e2 * Precision::Angular()
            bool radiusNormalToAxis = false;
            if (baseVec.Magnitude() > 1e2 * Precision::Angular())
                radiusNormalToAxis = baseVec.IsNormal(conSurface.Axis().Direction(), 1e2 * Precision::Angular());
            else if (topVec.Magnitude() > 1e2 * Precision::Angular())
                radiusNormalToAxis = topVec.IsNormal(conSurface.Axis().Direction(), 1e2 * Precision::Angular());
            if (radiusNormalToAxis)
            {
                // height is computed by finding projection of surface vector along the cone axis
                // cone height can be found using vectors:
                // e.g. u-periodic & alpha > 0:
                // auto center = gp_Vec(surface.Value(uFirst, vLast).XYZ()) + gp_Vec((0.5 * baseVec).XYZ());
                // auto apex = gp_Vec(surface.Value(uFirst, vFirst).XYZ());
                // auto height = (apex - center).Magnitude();
                // OC = OA + AC & CA = OA - OC (O = origin, C = base-center, A = apex, B = base-point)
                if (!semiArc)
                {
                    // cosine rule: computing radius from base vector
                    // due to cone symmetry, the same arc angle is applicable for base radius & top radius
                    // r = vec / (sqrt(2 * (1 - cos(theta))))
                    auto num = baseVec.Magnitude();
                    auto den = std::sqrt(2 * (1 - std::cos(arcAngle)));
                    auto rad = dimensions->radius = static_cast<float>(num / den);
                    num = topVec.Magnitude();
                    auto topRad = static_cast<float>(num / den);
                    auto hgt = fabs(surfVec.Dot(conSurface.Axis().Direction()));
                    dimensions->height = static_cast<float>(hgt);
                    // tan(alpha) = (base-rad - top-rad) / frustum-height
                    // if frustum is not present, top-rad = 0, frustum-height = height
                    // => tan(alpha) = base-rad/height
                    dimensions->angle = static_cast<float>(RAD2DEG(std::atan2(fabs(topRad - rad), hgt)));
                }
                else
                {
                    // normally CAD parts contain semi arc cones (base-angle = 180) and hence
                    // baseVec is along the diameter
                    auto topRad = 0.5 * topVec.Magnitude();
                    auto rad = 0.5 * baseVec.Magnitude();
                    auto hgt = fabs(surfVec.Dot(conSurface.Axis().Direction()));
                    dimensions->radius = rad > Precision::Confusion() ? static_cast<float>(rad) : static_cast<float>(topRad);
                    dimensions->height = static_cast<float>(hgt);
                    // tan(alpha) = (base-rad - top-rad) / frustum-height
                    // if frustum is not present, top-rad = 0, frustum-height = height
                    // => tan(alpha) = base-rad/height
                    dimensions->angle = static_cast<float>(RAD2DEG(std::atan2(fabs(topRad - rad), hgt)));
                }
                break;
            }
        }
    }
}


bool OCC_Computation::IdentifyUniqueDegenerateQuadrics( const BRepAdaptor_Surface& brepSurface,
                                                        std::vector<std::tuple<float, gp_Ax1, gp_Pnt>>& uniqueShapes)
{
    // identify unique holes: 2 or more faces could be a part of the same hole
    // hence radius, axis, center are checked to verify existing holes for a given face type
    // shapes are unique if any of the 3 conditions are met:
    // radii are different, axis/location are different, centers are different
    bool uniqueShape = true;

    gp_Pnt center = OCC_Computation::ComputeCenter(brepSurface);
    for (const auto& [rad, ax, ctr] : uniqueShapes)
    {
        if (GeomAbs_SurfaceType::GeomAbs_Cylinder == brepSurface.GetType())
        {
            gp_Cylinder cylSurface = brepSurface.Cylinder();
            uniqueShape &=
                (fabs(static_cast<float>(cylSurface.Radius()) - rad) > Precision::Confusion()) ||
                (!cylSurface.Axis().IsCoaxial(ax, Precision::Angular(), Precision::Confusion())) ||
                (!center.IsEqual(ctr, Precision::Confusion()));
        }
        else if (GeomAbs_SurfaceType::GeomAbs_Cone == brepSurface.GetType())
        {
            gp_Cone conSurface = brepSurface.Cone();
            uniqueShape &=
                (fabs(static_cast<float>(conSurface.RefRadius()) - rad) > Precision::Confusion()) ||
                (!conSurface.Axis().IsCoaxial(ax, Precision::Angular(), Precision::Confusion())) ||
                (!center.IsEqual(ctr, Precision::Confusion()));
        }
        else
            uniqueShape &= false;

        if (!uniqueShape)
            break;
    }

    if (uniqueShape)
    {
        if (GeomAbs_SurfaceType::GeomAbs_Cylinder == brepSurface.GetType())
        {
            // std::make_tuple()
            uniqueShapes.emplace_back(static_cast<float>(brepSurface.Cylinder().Radius()), \
                brepSurface.Cylinder().Axis(), center);
        }
        else if (GeomAbs_SurfaceType::GeomAbs_Cone == brepSurface.GetType())
        {
            // std::make_tuple()
            uniqueShapes.emplace_back(static_cast<float>(brepSurface.Cone().RefRadius()), \
                brepSurface.Cone().Axis(), center);
        }
    }

    return uniqueShape;
}


gp_Dir OCC_Computation::ComputePlaneNormal(const TopoDS_Shape& topoFace)
{
    auto brepSurface = BRepAdaptor_Surface(TopoDS::Face(topoFace));
    if (GeomAbs_SurfaceType::GeomAbs_Plane != brepSurface.GetType())
        return gp_Dir();
    gp_Dir faceNormal = brepSurface.Plane().Axis().Direction();
    if (TopAbs_Orientation::TopAbs_REVERSED == topoFace.Orientation())
        faceNormal.Reverse();
    return faceNormal;
}


bool OCC_Computation::IsCurveNonPlanar(const gp_Pnt* points)
{
    bool nonPlanar = true;
    auto vec1 = gp_Vec(points[1].XYZ()) - gp_Vec(points[0].XYZ());
    auto vec2 = gp_Vec(points[2].XYZ()) - gp_Vec(points[0].XYZ());
    auto vec3 = gp_Vec(points[3].XYZ()) - gp_Vec(points[0].XYZ());
    if (vec3.DotCross(vec1, vec2) < Precision::Confusion())
        nonPlanar = false;
    return nonPlanar;
}


bool OCC_Computation::IsInterpolatedCurveHelix( const BRepAdaptor_Surface& surface,
                                                const BRepAdaptor_Curve& curve,
                                                const TopoDS_Face& face,
                                                const TopoDS_Edge& edge,
                                                int16_t faceID, int16_t edgeID,
                                                HelixDimensions* dimensions)
{
    Handle(Geom_BSplineCurve) bspline = curve.BSpline();
    const TColgp_Array1OfPnt& poles = bspline->Poles();

    // 4 points are necessary to test coplanarity and helix is a space curve
    if (poles.Length() <= 3)
        return false;
    gp_Pnt points[4]{   bspline->StartPoint(),
                        poles(poles.Lower() + 1),
                        poles(poles.Upper() - 1),
                        bspline->EndPoint() };
    if (!OCC_Computation::IsCurveNonPlanar(points))
        return false;

    double low = 0.0, high = 0.0, surfaceRadius = 0.0;
    gp_Ax1 surfaceAxis;
    // 3 is order of derivative required
    // 1 -> tangent, 2 -> curvature/normal/binormal, 3 -> torsion
    GeomLProp_CLProps curveProps(BRep_Tool::Curve(edge, low, high), 3, Precision::Confusion());
    const TColStd_Array1OfReal& knots = bspline->Knots();
    auto meanRadius = 0.0f, meanPitch = 0.0f; auto count = 0;
    for (auto itr = knots.Lower(); itr <= knots.Upper(); ++itr)
    {
        // m = #knots, n = #poles, p = degree, k = knot-multiplicity
        // avoid the first & last points because the geometrical entities
        // behave differently (@k = p + 1)
        if (bspline->Multiplicity(itr) <= bspline->Degree())
        {
            curveProps.SetParameter(knots(itr));
            gp_Dir tangent, normal, binormal;
            curveProps.Tangent(tangent);
            curveProps.Normal(normal);
            binormal = tangent.Crossed(normal);

            // parent surface is the same irrespective of the iteration, compute only once
            if (!count++)
            {
                if (GeomAbs_Cone == surface.GetType())
                {
                    surfaceAxis = surface.Cone().Axis();
                    surfaceRadius = surface.Cone().RefRadius();
                }
                else
                {
                    surfaceAxis = surface.Cylinder().Axis();
                    surfaceRadius = surface.Cylinder().Radius();
                }
            }
            // tilt of helix is calculated as an acute angle wrt to the axis of parent surface
            // helix curvature vector is along the normal direction which is perpendicular to
            // the tilt, and accordingly, radius is computed horizontally and pitch vertically
            auto tilt = surfaceAxis.Direction().Angle(binormal);
            tilt = tilt > M_PI_2 ? M_PI - tilt : tilt;
            auto radius = (1 / curveProps.Curvature()) * std::cos(tilt);
            auto pitch = 2 * M_PI * radius * std::tan(tilt);
            meanRadius += radius;
            meanPitch += pitch;
#if 0
            const auto& d1 = curveProps.D1();
            const auto& d2 = curveProps.D2();
            const auto& d3 = curveProps.D3();
            // Lancret's theorem: curvature/torsion is constant for helix
            // torsion is numerically unstable, hence not computed
            auto torsion = d3.DotCross(d1, d2) / d1.CrossSquareMagnitude(d2);
            // curvature is identical to the previous calculation
            auto curvature = d1.CrossMagnitude(d2) / std::pow(d1.Magnitude(), 3);
            std::cout << "k: " << curvature << ", t: " << torsion << std::endl;
#endif
            if (Utility::DEBUG_LEVEL >= 3)
                std::cout << "Edge: " << edgeID << "k = " << knots(itr) << \
                ", a = " << radius << ", b = " << pitch << std::endl;
        }
    }

    meanRadius /= count;
    meanPitch /= count;
    // if the curve radius is within 10% tolerance of parent surface radius and
    // the curve has a positive pitch, the interpolated curve is a helix
    // checking the pitch could be redundant because test for coplanarity is
    // already done
    if ((fabs(surfaceRadius - meanRadius) / meanRadius) < 0.1f && meanPitch > Precision::Confusion())
    {
        dimensions->radius = meanRadius;
        dimensions->pitch = meanPitch;
        dimensions->parentRadius = surfaceRadius;
        dimensions->parentFace = faceID;
        dimensions->parentAxis = surfaceAxis;
        return true;
    }
    else
        return false;
}


void OCC_Computation::SortHelicesByAxes(const std::vector<HelixDimensions>& helices,
                                        std::vector<std::vector<HelixDimensions>>& sortedHelices)
{
    for (const auto& hitr : helices)
    {
        bool found = false;
        auto currentHelixAxis = hitr.parentAxis;
        for (auto& sitr : sortedHelices)
        {
            // parent surface radius is used because it is constant throughout the hole irrespective of
            // the helix curve radius (which may differ slightly)
            auto existingHelixAxis = sitr.at(0).parentAxis;
            auto existingHelixRadius = sitr.at(0).parentRadius;
            // 2 different helix curves can be collated (sorted) into different containers if they are spaced apart
            // by a distance greater than their parent surface radius (more likely) or if their axes directions are
            // different (less likely)
            if (currentHelixAxis.IsCoaxial(existingHelixAxis, 1e2 * Precision::Angular(), existingHelixRadius))
            {
                found = true;
                sitr.emplace_back(hitr.radius, hitr.pitch, hitr.parentRadius, hitr.parentFace, hitr.parentAxis);
                break;
            }
        }
        if (!found)
        {
            std::vector<HelixDimensions> temp{ HelixDimensions(
                hitr.radius, hitr.pitch, hitr.parentRadius, hitr.parentFace, hitr.parentAxis) };
            sortedHelices.emplace_back(std::move(temp));
        }
    }
}


void OCC_Computation::ComputeHelixModeDimensions(   const std::vector<HelixDimensions>& helixDimensions,
                                                    OCC_Computation::HelixDimensions* modeDimensions)
{
    // mode is computed because different helix curves within the same hole can have
    // slighly different pitch and/or radius and the one's which occur maximum times
    // should be chosen
    std::vector<float> helixParameters;
    helixParameters.reserve(helixDimensions.size());
    for (const auto& itr : helixDimensions)
        helixParameters.emplace_back(itr.radius);

    float modeRadius = 0.0f, modePitch = 0.0f;
    int16_t idxRadius = 0, idxPitch = 0, freqRadius = 0, freqPitch = 0;
    modeRadius = Utility::ComputeMode(helixParameters, idxRadius, freqRadius);

    // same vector (helixParameters) is used for computing mode of radius and pitch
    // reassign data from radius to pitch
    auto itrDim = helixDimensions.begin();
    for (auto& itrParam : helixParameters)
        itrParam = itrDim++->pitch;
    modePitch = Utility::ComputeMode(helixParameters, idxPitch, freqPitch);

    if (idxRadius != idxPitch)
    {
        if (freqRadius > freqPitch)
            modePitch = helixDimensions.at(idxRadius).pitch;
#if 0
        // radius has more relevance than pitch
        if (freqPitch > freqRadius)
            modeRadius = helixDimensions.at(idxPitch).radius;
#endif
    }

    *modeDimensions = helixDimensions.at(idxRadius);
};

