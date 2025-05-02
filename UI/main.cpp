#include "MyMainWindow.h"

int main (int theNbArgs, char** theArgVec)
{
  QApplication aQApp (theNbArgs, theArgVec);

  QCoreApplication::setApplicationName ("OCCT Qt Viewer");
  QCoreApplication::setOrganizationName ("OpenCASCADE");
  QCoreApplication::setApplicationVersion (OCC_VERSION_STRING_EXT);

#ifdef __APPLE__
  // suppress Qt warning "QCocoaGLContext: Falling back to unshared context"
  bool isCoreProfile = true;
  QSurfaceFormat aGlFormat;
  aGlFormat.setDepthBufferSize   (24);
  aGlFormat.setStencilBufferSize (8);
  if (isCoreProfile) { aGlFormat.setVersion (4, 5); }
  aGlFormat.setProfile (isCoreProfile ? QSurfaceFormat::CoreProfile : QSurfaceFormat::CompatibilityProfile);
  QSurfaceFormat::setDefaultFormat (aGlFormat);
#endif

  MyMainWindow aMainWindow;
  aMainWindow.resize (aMainWindow.sizeHint());
  aMainWindow.show();
  return aQApp.exec();
}
