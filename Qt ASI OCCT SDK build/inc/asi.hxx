#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <Standard_Handle.hxx>
#include <TopoDS_Shape.hxx>

#include <asiAlgo_AAG.h>

#include <map>
#include <vector>

class OCC_Computation;

class AS_FeatureRecognition
{
public:
    enum class PartCategory { UNDEFINED, MILLING, TURNING, LASER_CUTTING, GRINDING, GEAR_MANUFACTURING };
    AS_FeatureRecognition() = default;
    AS_FeatureRecognition(const char* cadFilePath);
    ~AS_FeatureRecognition();

    inline std::shared_ptr<OCC_Computation> GetCADPart() const
    { return this->m_CADPart; }
    inline Handle(asiAlgo_AAG) GetAAG() const
    { return this->m_AAG; }
    inline std::shared_ptr<std::map<int, TopoDS_Edge>> GetConvexEdges() const
    { return this->m_ConvexEdges; }
    inline std::shared_ptr<std::map<int, TopoDS_Edge>> GetConcaveEdges() const
    { return this->m_ConcaveEdges; }
    inline std::shared_ptr<std::map<int, TopoDS_Edge>> GetSmoothEdges() const
    { return this->m_SmoothC1Edges; }

    bool Initialize();
    bool STEPReader();
    bool IGESReader();

    enum class AS_FeatureRecognition::PartCategory IdentifyPartCategory();
    bool RunFeatureRecognition( const std::string& filePath,
                                int16_t iteration   );
    bool IdentifySlots(std::vector<std::array<float, 3>>& slotDimensions);
    bool IdentifyDrilledHoles(std::vector<std::pair<float, float>>& holeDimensions);
    bool IdentifyHandles(std::vector<std::pair<float, float>>& handleDimensions);

protected:
    bool GenerateAAG();
    bool ClassifyEdges(std::map<int, TopoDS_Edge>* undefinedEdges = nullptr);
    double ComputeDihedralAngle(const TopoDS_Face& face1,
                                const TopoDS_Face& face2,
                                bool allowC1Edges,
                                double angularTolerance,
                                std::map<int, TopoDS_Edge>* convexEdges = nullptr,
                                std::map<int, TopoDS_Edge>* concaveEdges = nullptr,
                                std::map<int, TopoDS_Edge>* smoothC1Edges = nullptr,
                                std::map<int, TopoDS_Edge>* undefinedEdges = nullptr);
    bool FindHoleFaces();
    bool ConvertToCanonicalGeometry();

private:
    std::shared_ptr<OCC_Computation> m_CADPart{ nullptr };

    Handle(asiAlgo_AAG) m_AAG;
    std::shared_ptr<std::map<int, TopoDS_Edge>> m_ConvexEdges{ nullptr };
    std::shared_ptr<std::map<int, TopoDS_Edge>> m_ConcaveEdges{ nullptr };
    std::shared_ptr<std::map<int, TopoDS_Edge>> m_SmoothC1Edges{ nullptr };
    std::shared_ptr<TColStd_PackedMapOfInteger> m_HoleFaces{ nullptr };
};
