#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

#include <Bnd_Box.hxx>
#include <Bnd_OBB.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <OSD_Timer.hxx>
#include <Standard_Version.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

#include <map>
#include <string>
#include <vector>

#define OCC 0


struct ConeCylinderDimensions
{
    float radius{ 0.0f }, height{ 0.0f }, angle{ 0.0f };
    gp_Pnt center;
};


// Handle(ClassName) = opencascade::handle<ClassName>
// to use handles with ClassName, that class must inherit from the Standard_Transient class
class OCC_Computation
{
public:
    OCC_Computation() = default;
    OCC_Computation(const char* cadFilePath,
                    TopoDS_Shape* shape = nullptr);
    ~OCC_Computation();

    int ReadFile(const char* cadFilePath = nullptr);
    int WriteFile(bool stepFormat = true);
    bool Compute(bool usingOCC = true);

    const std::string& GetFilePath() const
    { return this->m_FilePath; }
    inline std::shared_ptr<TopoDS_Shape> GetShape() const
    { return this->m_Shape; }
    inline void SetShape(std::shared_ptr<TopoDS_Shape>& shape)
    { this->m_Shape = shape; }
    bool GetTopoGeomMapsGenerated() const
    { return this->m_TopoGeomMapsGenerated; }
    int16_t GetGenus() const
    { return this->m_Genus; }
    int16_t GetHoles() const
    { return this->m_Holes; }
    float GetTolerance() const
    { return this->m_Tolerance; }

//protected:
    bool GenerateTopologyGeometryMaps();
    void ResetTopologyGeometryMaps();       // required only for IGES files
    bool QueryTopologyMaps( TopTools_IndexedMapOfShape& vertexMap,
                            TopTools_IndexedMapOfShape& edgeMap,
                            TopTools_IndexedMapOfShape& faceMap,
                            TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap,
                            TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap  );
    bool QueryGeometryMaps( std::map<int16_t, GeomAbs_CurveType>& edgeTypeMap,
                            std::map<int16_t, GeomAbs_SurfaceType>& faceTypeMap  );
    bool ConvertIGESToSolid();
    // 5 mm
    gp_XYZ ComputeBoundingBox(  float offset = 0.05f,
                                bool igesMode = false);
    void ComputePartGeometry(bool usingOCC = true);
    bool ComputePartGeometryMaps();
    std::tuple<int16_t, int16_t> ComputePartTopology();
    gp_Pnt ComputeCenter(const BRepAdaptor_Surface& surface);
    // Cylinders & Cones are degenerate quadrics
    // General Quadric Equation: Ax2 + By2 + Cz2 + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0
    // Cone: D = E = F = G = H = I = J = 0
    // Cylinder: A/B/C = 0, D = E = F = G = H = I = 0
    void ComputeDegenerateQuadricDimensions(const BRepAdaptor_Surface& surface,
                                            ConeCylinderDimensions* dimensions);
    bool IdentifyUniqueDegenerateQuadrics(  const BRepAdaptor_Surface& brepSurface,
                                            std::vector<std::tuple<float, gp_Ax1, gp_Pnt>>& uniqueShapes);
    gp_Dir ComputePlaneNormal(const TopoDS_Shape& topoFace);

private:
    std::string m_FilePath;
    std::shared_ptr<TopoDS_Shape> m_Shape{ nullptr };
    int16_t m_Genus{ 0 }, m_Holes{ 0 };
    bool m_TopoGeomMapsGenerated{ false };
    float m_Tolerance{ 0.0f };
    std::shared_ptr<Bnd_OBB> m_OrientedBoundingBox{ nullptr };
    std::shared_ptr<Bnd_Box> m_AxisAlignedBoundingBox{ nullptr };
    std::shared_ptr<TopTools_IndexedMapOfShape> m_VertexMap{ nullptr };
    std::shared_ptr<TopTools_IndexedMapOfShape> m_EdgeMap{ nullptr };
    std::shared_ptr<TopTools_IndexedMapOfShape> m_FaceMap{ nullptr };
    std::shared_ptr<TopTools_IndexedDataMapOfShapeListOfShape> m_VertexEdgeMap{ nullptr };
    std::shared_ptr<TopTools_IndexedDataMapOfShapeListOfShape> m_EdgeFaceMap{ nullptr };
    std::shared_ptr<std::map<int16_t, GeomAbs_CurveType>> m_EdgeTypeMap{ nullptr };
    std::shared_ptr<std::map<int16_t, GeomAbs_SurfaceType>> m_FaceTypeMap{ nullptr };
    std::shared_ptr<std::map<GeomAbs_CurveType, int16_t>> m_CurveGeometryMap{ nullptr };
    std::shared_ptr<std::map<GeomAbs_SurfaceType, int16_t>> m_SurfaceGeometryMap{ nullptr };
};
