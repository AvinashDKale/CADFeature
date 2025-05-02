// 3rd-party: OCC
#include <AIS_InteractiveContext.hxx>
#include <AIS_ViewController.hxx>
#include <Standard_Handle.hxx>
#include <Standard_Version.hxx>
#include <V3d_View.hxx>
// 3rd-party: QT
#include <QFileDialog>
#include <QListView>
#include <QOpenGLWidget>
#include <QPushButton>
#include <QStringList>
#include <QTreeView>


class AIS_ViewCube;
class OpenGl_Context;


//! Auxiliary wrapper to avoid OpenGL macros collisions between Qt and OCCT headers.
class OcctGlTools
{
public:
    //! Return GL context.
    static Handle(OpenGl_Context) GetGlContext(const Handle(V3d_View)& theView);
};


//! OCCT 3D View.
class OcctQtViewer : public QOpenGLWidget, public AIS_ViewController
{
    Q_OBJECT
public:

    //! Main constructor.
    OcctQtViewer(QWidget* theParent = nullptr);

    //! Destructor.
    virtual ~OcctQtViewer();

    //! Return Viewer.
    const Handle(V3d_Viewer)& Viewer() const { return myViewer; }

    //! Return View.
    const Handle(V3d_View)& View() const { return myView; }

    //! Return AIS context.
    const Handle(AIS_InteractiveContext)& Context() const { return myContext; }

    //! Return OpenGL info.
    const QString& getGlInfo() const { return myGlInfo; }

    //! Minial widget size.
    virtual QSize minimumSizeHint() const override { return QSize(200, 200); }

    //! Default widget size.
    virtual QSize sizeHint() const override { return QSize(1400, 800)/*QSize(720, 480)*/; }

#if (OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 8 && OCC_VERSION_MAINTENANCE >= 0)
    //! Handle subview focus change.
    virtual void OnSubviewChanged(  const Handle(AIS_InteractiveContext)&,
                                    const Handle(V3d_View)&,
                                    const Handle(V3d_View)& theNewView  ) override;
#endif

    std::shared_ptr<std::vector<std::string>> getFilePaths() const { return this->myFilePaths; }
    void setFilePaths(const std::vector<std::string>& filePaths) { *this->myFilePaths = filePaths; }

protected: // OpenGL events

    virtual void initializeGL() override;
    virtual void paintGL() override;
    //virtual void resizeGL(int , int ) override;

protected: // user input events

    virtual void closeEvent(QCloseEvent* theEvent) override;
    virtual void keyPressEvent(QKeyEvent* theEvent) override;
    virtual void mousePressEvent(QMouseEvent* theEvent) override;
    virtual void mouseReleaseEvent(QMouseEvent* theEvent) override;
    virtual void mouseMoveEvent(QMouseEvent* theEvent) override;
    virtual void wheelEvent(QWheelEvent* theEvent) override;

private:

    //! Dump OpenGL info.
    void dumpGlInfo(bool theIsBasic, bool theToPrint);

    //! Request widget paintGL() event.
    void updateView();

    //! Handle view redraw.
    virtual void handleViewRedraw(  const Handle(AIS_InteractiveContext)& theCtx,
                                    const Handle(V3d_View)& theView ) override;

    Handle(V3d_Viewer)             myViewer;
    Handle(V3d_View)               myView;
    Handle(V3d_View)               myFocusView;
    Handle(AIS_InteractiveContext) myContext;
    Handle(AIS_ViewCube)           myViewCube;

    QString myGlInfo;
    bool myIsCoreProfile;
    std::shared_ptr<std::vector<std::string>> myFilePaths{ nullptr };
};


class FileDialog : public QFileDialog
{
    Q_OBJECT
private:
    QListView* myListView;
    QTreeView* myTreeView;
    QPushButton* myBtnOpen;
    QStringList mySelectedFiles;

public slots:
    void chooseClicked();
public:
    FileDialog();
    QStringList selectedFiles();
    bool eventFilter(QObject* watched, QEvent* event);
};
