set(APPLICATION brdfExplorer)
project(${APPLICATION})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

set(HEADERS
	*.h)
set(SOURCES
	Quad.cpp
	Sphere.cpp
	BRDFAnalytic.cpp
	BRDFBase.cpp
	BRDFImageSlice.cpp
	BRDFMeasuredAniso.cpp
	BRDFMeasuredMERL.cpp
	ColorVarWidget.cpp
	FloatVarWidget.cpp
	DGLFrameBuffer.cpp
	DGLShader.cpp
	IBLWidget.cpp
	IBLWindow.cpp
	ImageSliceWidget.cpp
	ImageSliceWindow.cpp
	LitSphereWindow.cpp
	main.cpp
	glerror.cpp
	MainWindow.cpp
	ViewerWindow.cpp
	ParameterGroupWidget.cpp
	ParameterWindow.cpp
	SharedContextGLWidget.cpp
	ShowingDockWidget.cpp
	PlotCartesianWindow.cpp
	PlotCartesianWidget.cpp
	PlotPolarWidget.cpp
	Plot3DWidget.cpp
	LitSphereWidget.cpp
	SimpleModel.cpp
	Paths.cpp
	ptex/PtexReader.cpp
	ptex/PtexUtils.cpp
	ptex/PtexCache.cpp
	ptex/PtexHalf.cpp
)


find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
QT   += widgets opengl
DEFINES += PTEX_STATIC NOMINMAX

brdfs.path = $$DEST/share/brdf/brdfs
brdfs.files = ../brdfs/*

data.path = $$DEST/share/brdf/data
data.files = ../data/*

images.path = $$DEST/share/brdf/images
images.files = ../images/*

probes.path = $$DEST/share/brdf/probes
probes.files = ../probes/*

shaderTemplates.path = $$DEST/share/brdf/shaderTemplates
shaderTemplates.files = ../shaderTemplates/*

pkgconfig.path = $$DEST/$$LIBDIR/pkgconfig
pkgconfig.files = brdf.pc

INSTALLS = target brdfs data images probes shaderTemplates pkgconfig

