#pragma once

#include "graph.hpp"

//Parameterizing boundary of the mesh onto a disc
//Mode 0: Uniform Parameterization
//Mode 1: Chord Length Parameterization
void getBoundaryDiscParameterization(const matXd &Vb, matXd& Tb, double radius, int mode){ 
    double totalWeight = 0.0;
    vector<double> weights; 
    if(!mode){
        int n = Vb.rows(); 
        for(unsigned int i = 0; i < n; i++){
            weights.push_back(1.0 / n);
        }
    }else{
        for(unsigned int i = 0; i < Vb.rows(); i++){
            vec3d v1 = Vb.row(i);
            vec3d v2 = Vb.row((i + 1) % Vb.rows());
            weights.push_back((v2 - v1).norm());
            totalWeight += weights[i];
        }
        for(unsigned int i = 0; i < Vb.rows(); i++){
            weights[i] /= totalWeight;
        }
    }
    vector<vec2d> points; 
    double currentTheta = 0.0f; 
    unsigned int i = 0;
    while(i < Vb.rows()){
        points.push_back(vec2d(radius * cos(currentTheta), radius * sin(currentTheta)));
        currentTheta += 2 * M_PI * weights[i];
        i++;
    }
    Tb.resize(Vb.rows(), 2);
    for(unsigned int i = 0; i < Vb.rows(); i++){
        Tb.row(i) = points[i];
    }
}

void get3DParameterVertices(const matXd &Ti, matXd &To, vec3d centre, vec3d x, vec3d y){
    vector<vec3d> points;
    for(unsigned int i = 0; i < Ti.rows(); i++){
        vec3d point = centre + Ti(i, 0) * x + Ti(i, 1) * y;
        points.push_back(point);
    }
    To.resize(Ti.rows(), 3);
    for(unsigned int i = 0; i < Ti.rows(); i++){
        To.row(i) = points[i];
    }
}

tuple<vec3d, vec3d, vec3d, double> getDiscParameters(const matXd &V){
    vec3d centre(0,0,0);
    for(unsigned int i = 0; i < V.rows(); i++){
        centre += V.row(i);
    }
    centre /= V.rows();
    vec3d v1 = V.row(0).transpose() - centre;
    vec3d v2 = V.row(1).transpose() - centre;
    vec3d v3 = V.row(2).transpose() - centre;
    vec3d normal = v1.cross(v2);
    if(normal.dot(v3) < 0){
        normal *= -1;
    }
    vec3d x_axis = v1;
    vec3d y_axis = normal.cross(x_axis);
    double radius = max(v1.norm(), max(v2.norm(), v3.norm()));
    return make_tuple(centre, x_axis.normalized(), y_axis.normalized(), radius);
}