cmake_minimum_required(VERSION 3.20 FATAL_ERROR) 
set(CMAKE_VERBOSE_MAKEFILE ON)

project(AnalysisSitusApp C CXX)
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
find_package (OpenGL REQUIRED)

# Qt5
set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_AUTORCC ON)
set_property (GLOBAL PROPERTY AUTOMOC_SOURCE_GROUP "Generated Files/Moc")
set_property (GLOBAL PROPERTY AUTORCC_SOURCE_GROUP "Generated Files/Resources")
find_package (Qt5 COMPONENTS Widgets REQUIRED Core Gui Widgets OpenGL)

# glfw
find_package(glfw3 REQUIRED)

# Analysis-Situs
set(AS_SDK ${AS_PATH}/SDK)
include_directories(${AS_SDK}/include)
link_directories(${AS_SDK}/lib)
set(ASI_LIBS
	asiActiveData
	asiAlgo
	asiAsm
	asiData
	asiEngine
)

# project
include_directories(src inc)
set(SRC
	src/asi.cxx
	src/occ.cxx
	src/MyMainWindow.cpp
	src/OcctQtViewer.cpp
	src/OcctGlTools.cpp
	src/common.cxx
	src/FileDialog.cpp)
set(HDR
	inc/asi.hxx
	inc/occ.hxx
	inc/MyMainWindow.h
    inc/OcctQtViewer.h
    inc/OcctGlTools.h
	inc/common.hxx
	inc/FileDialog.h)

add_executable(${PROJECT_NAME} main.cpp ${SRC} ${HDR})
target_link_libraries(${PROJECT_NAME}
        ${OCCT_LIBS}
		${ASI_LIBS}
        Qt5::Core Qt5::Gui Qt5::Widgets Qt5::OpenGL
        glfw
        ${OPENGL_LIBRARIES}
)
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

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
# set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_ENVIRONMENT "\
    # PATH=%PATH%;${3RDPARTY_freetype_DLL_DIR};\
	# $<IF:$<CONFIG:Debug>,\
	# ${3RDPARTY_OCCT_DLL_DIR_DEBUG};${AS_SDK}/bind,\
	# ${3RDPARTY_OCCT_DLL_DIR};${AS_SDK}/bin>")
