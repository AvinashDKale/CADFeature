cmake_minimum_required(VERSION 3.20 FATAL_ERROR) 
set(CMAKE_VERBOSE_MAKEFILE ON)

project(CostEstimationTool C CXX)
set(CMAKE_CXX_STANDARD 17)

# 3rd party paths
set(AS_PATH "D:/CADFeatureRec/AnalysisSitus - Copy")
set(3RDPARTY_DIR ${AS_PATH}/analysissitus-libpack-msvc2019)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/asitus_3rdparty_macros.cmake)

# Eigen
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/asitus_Eigen.cmake)
include_directories(${3RDPARTY_EIGEN_DIR})

# OpenCASCADE
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/asitus_OCCT.cmake)
include (${CMAKE_SOURCE_DIR}/cmake/asitus_OCCT_3rdparty.cmake)
include_directories(${3RDPARTY_OCCT_INCLUDE_DIR})
link_directories(${3RDPARTY_OCCT_LIBRARY_DIR})
# extract all OCC libraries
set(OCCT_LIBS "")
file(GLOB LIB_PATHS ${3RDPARTY_OCCT_LIBRARY_DIR}/*)
foreach(LIB ${LIB_PATHS})
	# filename without file extension
	get_filename_component(LIB_NAME "${LIB}" NAME_WE)
	#message(STATUS ${LIB_NAME})
	list(APPEND OCCT_LIBS ${LIB_NAME})
endforeach()

# OpenGL
#set(OPENGL_DEBUG ${AS_PATH}/opengl-libpack/libd)
#set(OPENGL_RELEASE ${AS_PATH}/opengl-libpack/lib)
find_package (OpenGL REQUIRED)

# Qt
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set_property (GLOBAL PROPERTY AUTOMOC_SOURCE_GROUP "Generated Files/Moc")
set_property (GLOBAL PROPERTY AUTORCC_SOURCE_GROUP "Generated Files/Resources")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/asitus_Qt.cmake)
set(QT_LIBS
	Qt5::Core
	Qt5::Gui
	Qt5::Widgets
	Qt5::OpenGL)

# Analysis-Situs
set(AS_SDK ${AS_PATH}/SDK)
include_directories(${AS_SDK}/include)
link_directories(${AS_SDK}/lib)
set(ASI_LIBS
	asiActiveData
	asiAlgo
	asiAsm
	asiData
	asiEngine)

# project
include_directories(source include)
set(SRC
	source/asi.cxx
	source/common.cxx
	source/cost.cxx
	source/occ.cxx
	source/qt.cxx
	source/qt_main_window.cxx)
set(HDR
	include/asi.hxx
	include/common.hxx
	include/cost.hxx
	include/occ.hxx
	include/qt.hxx
	include/qt_main_window.hxx)

add_executable(${PROJECT_NAME} main.cpp ${SRC} ${HDR})
target_link_libraries(${PROJECT_NAME} ${OCCT_LIBS} ${QT_LIBS} ${ASI_LIBS})
target_link_libraries(${PROJECT_NAME} debug ${OPENGL_DEBUG}/GlU32.Lib ${OPENGL_DEBUG}/OpenGL32.Lib)
target_link_libraries(${PROJECT_NAME} optimized ${OPENGL_RELEASE}/GlU32.Lib ${OPENGL_RELEASE}/OpenGL32.Lib)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

# copy .dll/.lib/.pdb files of 3rd party pre-compiled library locations to AnalysisSitusApp.exe directory
# Analysis-Situs
add_custom_command(	TARGET ${PROJECT_NAME} POST_BUILD
					COMMAND ${CMAKE_COMMAND} -E copy_directory
					$<IF:$<CONFIG:Debug>,"${AS_SDK}/bind","${AS_SDK}/bin">
					$<TARGET_FILE_DIR:${PROJECT_NAME}>)
# OpenCASCADE
add_custom_command(	TARGET ${PROJECT_NAME} POST_BUILD
					COMMAND ${CMAKE_COMMAND} -E copy_directory
					$<IF:$<CONFIG:Debug>,"${3RDPARTY_OCCT_DLL_DIR_DEBUG}","${3RDPARTY_OCCT_DLL_DIR}">
					$<TARGET_FILE_DIR:${PROJECT_NAME}>)
# FreeType
add_custom_command(	TARGET ${PROJECT_NAME} POST_BUILD
					COMMAND ${CMAKE_COMMAND} -E copy_directory
					"${3RDPARTY_freetype_DLL_DIR}"
					$<TARGET_FILE_DIR:${PROJECT_NAME}>)
# Qt
add_custom_command(	TARGET ${PROJECT_NAME} POST_BUILD
					COMMAND ${CMAKE_COMMAND} -E copy_directory
					"${3RDPARTY_QT_DIR}/bin"
					$<TARGET_FILE_DIR:${PROJECT_NAME}>)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
set_target_properties (${PROJECT_NAME} PROPERTIES VS_DEBUGGER_ENVIRONMENT "QT_PLUGIN_PATH=${3RDPARTY_QT_DIR}/plugins;QT_QPA_PLATFORM_PLUGIN_PATH=${3RDPARTY_QT_DIR}/plugins/platforms")
