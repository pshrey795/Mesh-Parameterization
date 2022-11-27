#pragma once

//STL
#include <bits/stdc++.h>

//Libigl imports
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>

//OpenGL imports
#include <GLFW/glfw3.h>

using namespace std;
using namespace Eigen;

typedef MatrixXd matXd;
typedef MatrixXi matXi;
typedef VectorXd vecXd;
typedef VectorXi vecXi;
typedef Vector2d vec2d;
typedef Vector3d vec3d;

#define OBJECT_PATH "../assets/objects/"

//Global Variables
//Mesh
matXd V; matXi F;

//GUI
igl::opengl::glfw::Viewer viewer;
igl::opengl::glfw::imgui::ImGuiPlugin plugin;
igl::opengl::glfw::imgui::ImGuizmoWidget gizmo;
igl::opengl::glfw::imgui::ImGuiMenu menu;
