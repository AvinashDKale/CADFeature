#include <array>
#include <map>
#include <vector>
// 3rd-party: ASI
#include <asiAlgo_AAG.h>
// 3rd-party: OCC
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Standard_Handle.hxx>
#include <TopoDS_Shape.hxx>


class OCC_Computation;

class AS_FeatureRecognition
{
public:
    AS_FeatureRecognition() = default;
    explicit AS_FeatureRecognition( const char* cadFilePath,
                                    TopoDS_Shape* shape = nullptr  );
    virtual ~AS_FeatureRecognition();

    // prohibit lvalue copy/assignment since FR is specific to a CAD part
    AS_FeatureRecognition(const AS_FeatureRecognition&) = delete;
    AS_FeatureRecognition& operator=(const AS_FeatureRecognition&) = delete;
    void MoveSemantics(AS_FeatureRecognition& rhs) noexcept(true);
    AS_FeatureRecognition(AS_FeatureRecognition&& rhs) noexcept(true);
    AS_FeatureRecognition& operator=(AS_FeatureRecognition&& rhs) noexcept(true);

    inline std::shared_ptr<OCC_Computation> GetCADPart() const noexcept(true)
    { return this->m_CADPart; }
    inline Handle(asiAlgo_AAG) GetAAG() const noexcept(true)
    { return this->m_AAG; }
    inline std::shared_ptr<std::map<int16_t, TopoDS_Edge>> GetConvexEdges() const noexcept(true)
    { return this->m_ConvexEdges; }
    inline std::shared_ptr<std::map<int16_t, TopoDS_Edge>> GetConcaveEdges() const noexcept(true)
    { return this->m_ConcaveEdges; }
    inline std::shared_ptr<std::map<int16_t, TopoDS_Edge>> GetSmoothEdges() const noexcept(true)
    { return this->m_SmoothC1Edges; }
    inline std::shared_ptr<std::map<int16_t, TopoDS_Edge>> GetHoleEdges() const noexcept(true)
    { return this->m_HoleEdges; }

    // cad file readers
    bool STEPReader();
    bool IGESReader();

    // feature recognition APIs
    bool RunFeatureRecognition( const std::string& filePath,
                                std::vector<std::array<float, 3>>& slotDimensions,
                                std::vector<std::pair<float, float>>& drilledHoleDimensions,
                                std::vector<std::pair<float, float>>& threadedHoleDimensions,
                                std::vector<std::pair<float, float>>& handleDimensions,
                                std::array<float, 3>& boxDimensions );
    bool IdentifySlots(std::vector<std::array<float, 3>>& slotDimensions) const;
    bool IdentifyHoles( bool drilledHoles,
                        std::vector<std::pair<float, float>>& holeDimensions) const;
    bool IdentifyHandles(std::vector<std::pair<float, float>>& handleDimensions) const;

protected:
    enum class PartCategory { UNDEFINED, MILLING, TURNING, LASER_CUTTING, GRINDING, GEAR_MANUFACTURING };
    enum class AS_FeatureRecognition::PartCategory IdentifyPartCategory();
    bool Initialize();
    bool GenerateAAG();
    bool ClassifyEdges(std::map<int16_t, TopoDS_Edge>* undefinedEdges = nullptr);
    double ComputeDihedralAngle(const TopoDS_Face& face1,
                                const TopoDS_Face& face2,
                                bool allowC1Edges,
                                double angularTolerance,
                                std::map<int16_t, TopoDS_Edge>* convexEdges = nullptr,
                                std::map<int16_t, TopoDS_Edge>* concaveEdges = nullptr,
                                std::map<int16_t, TopoDS_Edge>* smoothC1Edges = nullptr,
                                std::map<int16_t, TopoDS_Edge>* undefinedEdges = nullptr) const;
    bool ComputeHoleFaces();
    bool ComputeHoleEdges();
    bool ComputeShaftFaces();
    bool ConvertToCanonicalGeometry();
    bool InvertShellOrientation();
    bool AdjacentHelicoidPresent(t_topoId faceID) const;
    void PrintFaceIDGeometry(bool featureCategory) const;

private:
    std::shared_ptr<OCC_Computation> m_CADPart{ nullptr };

    Handle(asiAlgo_AAG) m_AAG;
    std::shared_ptr<std::map<int16_t, TopoDS_Edge>> m_ConvexEdges{ nullptr };
    std::shared_ptr<std::map<int16_t, TopoDS_Edge>> m_ConcaveEdges{ nullptr };
    std::shared_ptr<std::map<int16_t, TopoDS_Edge>> m_SmoothC1Edges{ nullptr };
    std::shared_ptr<std::map<int16_t, TopoDS_Edge>> m_HoleEdges{ nullptr };
    std::shared_ptr<std::map<int16_t, std::pair<float, float>>> m_Threads{ nullptr };
    std::shared_ptr<TColStd_PackedMapOfInteger> m_HoleFaces{ nullptr };
    std::shared_ptr<TColStd_PackedMapOfInteger> m_ShaftFaces{ nullptr };
    // inline keyword on static variables obviates the need to define them in source file
    inline static int count{ 0 };
    bool m_ThreadsPresent{ false };
};
