#include <array>
#include <cassert>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>
// 3rd-party: OCC
#include <Bnd_Box.hxx>
#include <Bnd_OBB.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>
#include <OSD_Timer.hxx>
#include <Standard_Version.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>


// Handle(ClassName) = opencascade::handle<ClassName>
// to use handles with ClassName, that class must inherit from the Standard_Transient class
class OCC_Computation
{
public:
    struct ConeCylinderDimensions
    {
        float radius{ 0.0f }, height{ 0.0f }, angle{ 0.0f };
        gp_Pnt center;
    };
    struct HelixDimensions
    {
        float radius{ 0.0f }, pitch{ 0.0f }, parentRadius{ 0.0f };
        int16_t parentFace{ 0 };
        gp_Ax1 parentAxis;
        HelixDimensions() = default;
        HelixDimensions(float param1, float param2, float param3, int16_t param4, gp_Ax1 param5) : 
            radius(param1), pitch(param2), parentRadius(param3), parentFace(param4), parentAxis (param5) {}
        ~HelixDimensions() = default;
    };

    OCC_Computation(const char* cadFilePath = nullptr,
                    TopoDS_Shape* shape = nullptr);
    virtual ~OCC_Computation();

    // prohibit lvalue copy/assignment
    OCC_Computation(const OCC_Computation&) = delete;
    OCC_Computation& operator=(const OCC_Computation&) = delete;
    void MoveSemantics(OCC_Computation& rhs) noexcept(true);
    OCC_Computation(OCC_Computation&& rhs) noexcept(true);
    OCC_Computation& operator=(OCC_Computation&& rhs) noexcept(true);

    // cad file i/o
    // format: STEP = 0, IGES = 1, BREP = 2
    int ReadFile(const char* cadFilePath = nullptr);
    int WriteFile(int16_t format = 0);

    // cad computation
    bool Compute(   bool usingOCC = true,
                    gp_XYZ* boundingBox = nullptr   );
    // standard offset is 5 mm (for the assumed dimension range of parts)
    gp_XYZ ComputeBoundingBox(  float offset = 0.05f,
                                bool igesMode = false   );
    void ComputePartGeometry(bool usingOCC = true) const;
    std::pair<int16_t, int16_t> ComputePartTopology();

    const std::string& GetFilePath() const noexcept(true)
    { return this->m_FilePath; }
    inline std::shared_ptr<TopoDS_Shape> GetShape() const noexcept(true)
    { return this->m_Shape; }
    inline void SetShape(std::shared_ptr<TopoDS_Shape>& shape)
    { /*assert(!shape->IsNull());*/ this->m_Shape = shape; }
    constexpr bool GetTopoGeomMapsGenerated() const noexcept(true)
    { return this->m_TopoGeomMapsGenerated; }
    constexpr int16_t GetGenus() const noexcept(true)
    { return this->m_Genus; }
    constexpr int16_t GetHoles() const noexcept(true)
    { return this->m_Holes; }
    constexpr float GetTolerance() const noexcept(true)
    { return this->m_Tolerance; }

    // TBD : change the class design (HLD) to keep the below methods protected in both
    // OCC_Computation & AS_FeatureRecognition classes
//protected:
    bool ComputeCurveSurfaceGeometryMaps();
    bool GenerateTopologyGeometryMaps();
    void ResetTopologyGeometryMaps();       // required primarily for IGES files
    bool UpdateTopologyGeometryMaps();
    bool QueryTopologyMaps( TopTools_IndexedMapOfShape& vertexMap,
                            TopTools_IndexedMapOfShape& edgeMap,
                            TopTools_IndexedMapOfShape& faceMap,
                            TopTools_IndexedDataMapOfShapeListOfShape& vertexEdgeMap,
                            TopTools_IndexedDataMapOfShapeListOfShape& edgeFaceMap  ) const;
    bool QueryGeometryMaps( std::map<int16_t, GeomAbs_CurveType>& edgeTypeMap,
                            std::map<int16_t, GeomAbs_SurfaceType>& faceTypeMap  ) const;
    bool ConvertIGESToSolid();
    std::pair<int16_t, int16_t> CountProceduralGeometryEntities();
    std::pair<int16_t, int16_t> ConvertProceduralGeometryToBSplineGeometry();
    std::array<int16_t, 3> FixRedundantTopology();

    // Cylinders & Cones are degenerate quadrics
    // General Quadric Equation: Ax2 + By2 + Cz2 + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0
    // Cone: D = E = F = G = H = I = J = 0
    // Cylinder: A/B/C = 0, D = E = F = G = H = I = 0
    static const char* CurveGeometry(GeomAbs_CurveType curveType);
    static const char* SurfaceGeometry(GeomAbs_SurfaceType surfaceType);
    static gp_Pnt ComputeCenter(const BRepAdaptor_Surface& surface);
    static void ComputeDegenerateQuadricDimensions( const BRepAdaptor_Surface& surface,
                                                    ConeCylinderDimensions* dimensions  );
    static bool IdentifyUniqueDegenerateQuadrics(   const BRepAdaptor_Surface& brepSurface,
                                                    std::vector<std::tuple<float, gp_Ax1, gp_Pnt>>& uniqueShapes);
    static gp_Dir ComputePlaneNormal(const TopoDS_Shape& topoFace);
    static bool IsCurveNonPlanar(const gp_Pnt *points);
    static bool IsInterpolatedCurveHelix(   const BRepAdaptor_Surface& surface,
                                            const BRepAdaptor_Curve& curve,
                                            const TopoDS_Face& face,
                                            const TopoDS_Edge& edge,
                                            int16_t faceID, int16_t edgeID,
                                            HelixDimensions* dimensions );
    static void SortHelicesByAxes(  const std::vector<HelixDimensions>& helices,
                                    std::vector<std::vector<HelixDimensions>>& sortedHelices);
    static void ComputeHelixModeDimensions( const std::vector<HelixDimensions>& helixDimensions,
                                            HelixDimensions* modeDimensions);

private:
    std::string m_FilePath;
    mutable std::shared_ptr<TopoDS_Shape> m_Shape{ nullptr };
    // m_Genus is a through hole, m_Holes is a blind hole
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
