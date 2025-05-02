#include <common.hxx>
#include <occ.hxx>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
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
#include <Geom_ElementarySurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GProp_GProps.hxx>
#include <IGESControl_Reader.hxx> 
#include <IGESControl_Writer.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <ShapeAnalysis_ShapeContents.hxx>
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
#include <XSControl_Reader.hxx>
#include <XSControl_Writer.hxx>


OCC_Computation::OCC_Computation(   const char* cadFilePath,
                                    TopoDS_Shape* shape )
{
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
    //this->m_CurveGeometryMap = std::make_shared<std::map<GeomAbs_CurveType, int16_t>>();
    //this->m_SurfaceGeometryMap = std::make_shared<std::map<GeomAbs_SurfaceType, int16_t>>();
    //std::cout << "Constructor!" << std::endl;
}


OCC_Computation::~OCC_Computation()
{
    //std::cout << "Destructor!" << std::endl;
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
        //auto edgeID = ;
        //BRepAdaptor_Curve curve(TopoDS::Edge((*this->m_EdgeMap)(idx)));
        this->m_EdgeTypeMap->emplace(std::make_pair(idx, BRepAdaptor_Curve(TopoDS::Edge((*this->m_EdgeMap)(idx))).GetType()));
    }
    for (int16_t idx = 1; idx <= this->m_FaceMap->Extent(); ++idx)
    {
        //auto faceID = ;
        //BRepAdaptor_Surface surface(TopoDS::Face(faceID));
        this->m_FaceTypeMap->emplace(std::make_pair(idx, BRepAdaptor_Surface(TopoDS::Face((*this->m_FaceMap)(idx))).GetType()));
    }

    this->m_TopoGeomMapsGenerated = true;
    return true;
}


void OCC_Computation::ResetTopologyGeometryMaps()
{
    this->m_VertexMap->Clear();
    this->m_EdgeMap->Clear();
    this->m_FaceMap->Clear();
    this->m_VertexEdgeMap->Clear();
    this->m_EdgeFaceMap->Clear();
    this->m_EdgeTypeMap->clear();
    this->m_FaceTypeMap->clear();
    this->m_TopoGeomMapsGenerated = false;
}


bool OCC_Computation::QueryTopologyMaps(TopTools_IndexedMapOfShape& vertexMap,
                                        TopTools_IndexedMapOfShape& edgeMap,
                                        TopTools_IndexedMapOfShape& faceMap,
                                        TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap,
                                        TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap  )
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
                                        std::map<int16_t, GeomAbs_SurfaceType>& faceTypeMap)
{
    if (!this->m_VertexMap || this->m_VertexMap->IsEmpty()) return false;
    edgeTypeMap = *this->m_EdgeTypeMap;
    faceTypeMap = *this->m_FaceTypeMap;
    return true;
}


bool OCC_Computation::Compute(bool usingOCC)
{
    bool status = true;
    if (!this->m_Shape->IsNull())
    {
        std::cout << "\nCAD Computation:" << std::endl;
        //if (this->m_AxisAlignedBoundingBox->SquareExtent() < Precision::SquareConfusion())
        gp_XYZ boundingBox = this->ComputeBoundingBox();
        if (usingOCC && !this->m_TopoGeomMapsGenerated)
            this->GenerateTopologyGeometryMaps();
        this->ComputePartTopology();
        //this->ComputePartGeometryMaps();
        if (usingOCC)
            this->ComputePartGeometry();
    }
    else
        status = false;
    return status;
}


int OCC_Computation::ReadFile(const char* cadFilePath)
{
    // [&] captures this* and all local variables
    auto ReadOperation = [&](XSControl_Reader* reader) -> int
    {
        IFSelect_ReturnStatus status = reader->ReadFile(cadFilePath);
        if (IFSelect_RetDone != status)
            return -1;
        int numTranslations = reader->TransferRoots();
        if (!numTranslations)
            return -1;
        *this->m_Shape = reader->OneShape();
        if (this->m_Shape->IsNull())
            return -1;
        return 0;
    };

    if (!cadFilePath)
        cadFilePath = this->m_FilePath.c_str();
    int status = 0;
    bool igesFile = false;
    auto fileExtn = Utility::FileExtn(cadFilePath);
    if ("step" == fileExtn || "stp" == fileExtn)
    {
        STEPControl_Reader stpReader;
        status = ReadOperation(&stpReader);
    }
    else if ("iges" == fileExtn || "igs" == fileExtn)
    {
        IGESControl_Reader igsReader;
        status = ReadOperation(&igsReader);
        igesFile = true;
    }
    else
    {
        return -1;
    }

    if (status)
        return -1;
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

    if (igesFile)
    {
        status &= this->GenerateTopologyGeometryMaps() ? 0 : -1;
        status &= this->ConvertIGESToSolid() ? 0 : -1;
        // necessary because conversion of surfaces to solids changes vertex-edge & edge-face maps
        this->ResetTopologyGeometryMaps();
    }
    std::cout << "\n\nShapeType: " << TopAbs::ShapeTypeToString(this->m_Shape->ShapeType()) << std::endl;

    return status;
}


int OCC_Computation::WriteFile(bool stepFormat)
{
    int status = 0;
    if (this->m_FilePath.empty() || this->m_Shape->IsNull())
        return -1;

    auto filePath = this->m_FilePath.c_str();
    auto fileExtn = Utility::FileExtn(filePath);
    if ((stepFormat && ("step" == fileExtn || "stp" == fileExtn)) ||
        (!stepFormat && ("iges" == fileExtn || "igs" == fileExtn)))
        return -1;
    auto outputPath = Utility::OutputPath(filePath, stepFormat);

    if (stepFormat)
    {
        STEPControl_Writer stpWriter;
        IFSelect_ReturnStatus status = stpWriter.Transfer(*this->m_Shape, STEPControl_ManifoldSolidBrep);
        if (IFSelect_RetDone != status)
            status = stpWriter.Transfer(*this->m_Shape, STEPControl_BrepWithVoids);
        if (IFSelect_RetDone != status)
            return -1;
        status = stpWriter.Write(outputPath.c_str());
        if (IFSelect_RetDone != status)
            return -1;
    }
    else
    {
        IGESControl_Writer igsWriter;
        bool status = igsWriter.AddShape(*this->m_Shape);
        if (!status)
            return -1;
        status = igsWriter.Write(outputPath.c_str());
        if (!status)
            return -1;
    }

    return 0;
}


bool OCC_Computation::ConvertIGESToSolid()
{
    auto MinMaxAvgeTolerance = [](  const std::shared_ptr<TopTools_IndexedMapOfShape>& topoMap,
                                    int topoType, double& min, double& max, double& avge) -> void
    {
        for (int idx = 1; idx <= topoMap->Extent(); ++idx)
        {
            auto ithEntity = (*topoMap)(idx);
            auto ithTolerance = 0.0;
            if (!topoType)
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

    // return if the shape is already a solid
    if (TopAbs_SOLID == this->m_Shape->ShapeType())
       return false;

    double minVertexTolerance = DBL_MAX, minEdgeTolerance = DBL_MAX, avgeVertexTolerance = 0.0;
    double maxVertexTolerance = DBL_MIN, maxEdgeTolerance = DBL_MIN, avgeEdgeTolerance = 0.0;
    MinMaxAvgeTolerance(this->m_VertexMap, 0, minVertexTolerance, maxVertexTolerance, avgeVertexTolerance);
    MinMaxAvgeTolerance(this->m_EdgeMap, 1, minEdgeTolerance, maxEdgeTolerance, avgeEdgeTolerance);
    //auto maxFaceTolerance = BRep_Tool::MaxTolerance(*this->m_Shape, TopAbs_FACE);

    auto minTol = std::min(avgeVertexTolerance, avgeEdgeTolerance);
    auto maxTol = std::max(1e1 * maxVertexTolerance, 1e1 * maxEdgeTolerance);
    auto boundingBoxDiagonal = this->ComputeBoundingBox(0.0f, true).Modulus();
    // criteria can be modified
    this->m_Tolerance = std::min(maxTol, 1e-3 * boundingBoxDiagonal);

    BRepBuilderAPI_Sewing sewFaces(this->m_Tolerance);
    sewFaces.SetMinTolerance(minTol);
    sewFaces.SetMaxTolerance(maxTol);
    for (TopExp_Explorer explorer(*this->m_Shape, TopAbs_FACE); explorer.More(); explorer.Next())
        sewFaces.Add(explorer.Current());
    sewFaces.Perform();
    const TopoDS_Shape& sewedShape = sewFaces.SewedShape();
    // if the shape is not a shell of faces, tolerance needs to be modified
    if (TopAbs_SHELL != sewedShape.ShapeType())
        return false;

    BRepBuilderAPI_MakeSolid solid(TopoDS::Shell(sewedShape));
    this->m_Shape = std::make_shared<TopoDS_Shape>(solid.Shape());

   return true;
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
    static int count = 0;
    if (!igesMode || !count++)
    {
        // compute AABB only once when for IGES parts
        BRepBndLib::Add(*this->m_Shape, *this->m_AxisAlignedBoundingBox);
        ComputeAxesDim(*this->m_AxisAlignedBoundingBox, axesDims);
    }

    gp_XYZ offsettedBoundingBox(axesDims[0], axesDims[1], axesDims[2]);
    if (!igesMode)
    {
        // query AABB dimensions before enlarging the box
        ComputeAxesDim(*this->m_AxisAlignedBoundingBox, axesDims);

        // PCA is used for finding PMI (principal moment of inertia) axes
        BRepBndLib::AddOBB(*this->m_Shape, *this->m_OrientedBoundingBox, false, true, true);
        // OBB generates box center and half dimensions along each axes
        pcaDims[0] = 2.0f * this->m_OrientedBoundingBox->XHSize();
        pcaDims[1] = 2.0f * this->m_OrientedBoundingBox->YHSize();
        pcaDims[2] = 2.0f * this->m_OrientedBoundingBox->ZHSize();

        // offsetted bounding box
        this->m_AxisAlignedBoundingBox->Enlarge(offset);
        ComputeAxesDim(*this->m_AxisAlignedBoundingBox, offAxesDims);
        offsettedBoundingBox.SetX(offAxesDims[0]);
        offsettedBoundingBox.SetY(offAxesDims[1]);
        offsettedBoundingBox.SetZ(offAxesDims[2]);

        std::cout << "Axis-Aligned-Bounding-Box: (" << axesDims[0] << ", " << axesDims[1] << ", " << axesDims[2] << ")" << std::endl;
        std::cout << "Oriented-Bounding-Box: (" << pcaDims[0] << ", " << pcaDims[1] << ", " << pcaDims[2] << ")" << std::endl;
        //std::cout << "OBB == AABB: " << std::boolalpha << this->m_OrientedBoundingBox->IsAABox() << std::endl;
        std::cout << "Offsetted-Bounding-Box: (" << \
            offsettedBoundingBox.X() << ", " << offsettedBoundingBox.Y() << ", " << offsettedBoundingBox.Z() << ")" << std::endl;
    }

    return offsettedBoundingBox;
}


void OCC_Computation::ComputePartGeometry(bool usingOCC)
{
    auto CurveType = [](GeomAbs_CurveType type) -> const char*
    {
        switch (type)
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
    };

    auto SurfaceType = [](GeomAbs_SurfaceType type) -> const char*
    {
        switch (type)
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
    };

    if (this->m_VertexMap->IsEmpty())
        return;

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
        auto crvType = CurveType(curve.GetType());
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
        auto srfType = SurfaceType(surface.GetType());
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


bool OCC_Computation::ComputePartGeometryMaps()
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

    if (this->m_VertexMap->IsEmpty())
        return false;

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


std::tuple<int16_t, int16_t> OCC_Computation::ComputePartTopology()
{
    if (this->m_EdgeMap->IsEmpty())
        return { -1, -1 };

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
    std::cout << "\nCAD Topology:" << std::endl;
    topology.DumpExtent(std::cout);
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
    std::cout << "DEGENERECIES: " << numDegenerateEdges << ", MULTI-LOOPS = " << numMultiLoops << \
        ", HOLES <= " << this->m_Holes << ", GENUS = " << this->m_Genus << "\n\n";
    return { numDegenerateEdges, numMultiLoops };
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
    gp_Vec surfVec, baseVec;

    // surface vector on the cylinder surface is formed such that it is parallel to the axis
    // base vector is formed such that it is along the base radius
    for (int16_t i = 0; i < 2; ++i)
    {
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
                    // vLast = 0 (circular-base), vFirst < 0 (apex)
                    if (!i)
                        surfVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uFirst, vFirst)));
                    else
                        surfVec.Add(gp_Vec(surface.Value(uLast, vLast), surface.Value(uFirst, vFirst)));
                    baseVec.Add(gp_Vec(surface.Value(uFirst, vLast), surface.Value(uLast, vLast)));
                }
                else
                {
                    // vLast = 0 (circular-base), vFirst > 0 (apex) - to be implemented
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
                    // uLast = 0 (circular-base), uFirst < 0 (apex)
                    if (!i)
                        surfVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uFirst, vFirst)));
                    else
                        surfVec.Add(gp_Vec(surface.Value(uLast, vLast), surface.Value(uFirst, vFirst)));
                    baseVec.Add(gp_Vec(surface.Value(uLast, vFirst), surface.Value(uLast, vLast)));
                }
                else
                {
                    // uLast = 0 (circular-base), uFirst > 0 (apex) - to be implemented
                    ;
                }
            }
        }

        if (GeomAbs_SurfaceType::GeomAbs_Cylinder == surface.GetType())
        {
            gp_Cylinder cylSurface = surface.Cylinder();
            dimensions->center = this->ComputeCenter(surface);
            if (surfVec.IsParallel(cylSurface.Axis().Direction(), Precision::Angular()) &&
                (baseVec.IsNormal(cylSurface.XAxis().Direction(), Precision::Angular()) ||
                baseVec.IsNormal(cylSurface.YAxis().Direction(), Precision::Angular())))
            {
                dimensions->radius = static_cast<float>(0.5 * baseVec.Magnitude());
                dimensions->height = static_cast<float>(surfVec.Magnitude());
                break;
            }
        }
        else if (GeomAbs_SurfaceType::GeomAbs_Cone == surface.GetType())
        {
            gp_Cone conSurface = surface.Cone();
            if (baseVec.IsNormal(conSurface.XAxis().Direction(), Precision::Angular()) ||
                baseVec.IsNormal(conSurface.YAxis().Direction(), Precision::Angular()))
            {
                // height is computed by finding projection of surface vector along the cone axis
                // cone height can be found using vectors:
                // e.g. u-periodic & alpha > 0:
                // auto center = gp_Vec(surface.Value(uFirst, vLast).XYZ()) + gp_Vec((0.5 * baseVec).XYZ());
                // auto apex = gp_Vec(surface.Value(uFirst, vFirst).XYZ());
                // auto height = (apex - center).Magnitude();
                dimensions->radius = static_cast<float>(0.5 * baseVec.Magnitude());
                dimensions->height = fabs(surfVec.Dot(conSurface.Axis().Direction()));
                dimensions->angle = RAD2DEG(std::atan2f(dimensions->radius, dimensions->height));
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

    gp_Pnt center = this->ComputeCenter(brepSurface);
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
            uniqueShapes.emplace_back(std::make_tuple(static_cast<float>(brepSurface.Cylinder().Radius()), \
                brepSurface.Cylinder().Axis(), center));
        }
        else if (GeomAbs_SurfaceType::GeomAbs_Cone == brepSurface.GetType())
        {
            uniqueShapes.emplace_back(std::make_tuple(static_cast<float>(brepSurface.Cone().RefRadius()), \
                brepSurface.Cone().Axis(), center));
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
