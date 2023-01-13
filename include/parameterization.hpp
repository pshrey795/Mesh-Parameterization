#pragma once

#include "common.hpp"

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
    To.resize(Ti.rows(), 3);
    for(unsigned int i = 0; i < Ti.rows(); i++){
        To.row(i) = centre + Ti(i, 0) * x + Ti(i, 1) * y;
    }
}

tuple<vec3d, vec3d, vec3d, double> getDiscParameters(const matXd &V){
    vec3d centre(0.0,0.0,0.0);
    for(unsigned int i = 0; i < V.rows(); i++){
        centre += V.row(i);
    }
    centre /= V.rows();
    double radius = 0.0f;
    for(unsigned int i = 0; i < V.rows(); i++){
        radius = max(radius, (V.row(i).transpose() - centre).norm());
    }
    vec3d v1 = V.row(0).transpose() - centre;
    vec3d v2 = V.row(1).transpose() - centre;
    vec3d normal = v1.cross(v2);
    vec3d x_axis = v1;
    vec3d y_axis = normal.cross(x_axis);
    return make_tuple(centre, x_axis.normalized(), y_axis.normalized(), radius);
}

vector<int> getPinnedPoints(const matXd &V, const vecXi &U, int mode){
    if(mode){
        double maxDist = 0.0f;
        pair<int, int> farthestPair; 
        for(int i = 0; i < U.rows(); i++){
            for(int j = i + 1; j < U.rows(); j++){
                double currentDist = (V.row(U(i)) - V.row(U(j))).norm();
                if(currentDist > maxDist){
                    maxDist = currentDist;
                    farthestPair = make_pair(U(i), U(j));
                }
            }
        } 
        vector<int> v = {farthestPair.first, farthestPair.second};
        return v;
    }else{
        double maxDist = 0.0f;
        pair<int, int> farthestPair; 
        for(int i = 0; i < V.rows(); i++){
            for(int j = i + 1; j < V.rows(); j++){
                double currentDist = (V.row(i) - V.row(j)).norm();
                if(currentDist > maxDist){
                    maxDist = currentDist;
                    farthestPair = make_pair(i, j);
                }
            }
        } 
        vector<int> v = {farthestPair.first, farthestPair.second};
        return v;
    }
}

void fillMatrix(matXd &A, matXd &B, int vertexIdx, int faceIdx, int faceSz, int vertexSz, int pinnedSz, double Wr, double Wi, const vector<int> &pinnedPoints, const unordered_map<int, int> &pinnedPointsMap, const unordered_map<int, int> &effectiveIndexMap){
    int i = faceIdx;
    int j = vertexIdx;
    int n_dash = faceSz;
    int n = vertexSz;
    int p = pinnedSz;

    if(pinnedPointsMap.find(j) == pinnedPointsMap.end()){
        //Free vertex, so contributions will be added to LHS i.e. A
        //Entries (i,j'), (n' + i, j'), (i, n - p + j'), (n' + i, n - p + j') will be modified
        //where j' = Effective index of current vertex after subtracting the number of pinned points 
        //whose index is less than j
        j = effectiveIndexMap.at(j);
        A(i, j) = Wr; A(i + n_dash, j + n - p) = Wr;
        A(i, j + n - p) = -Wi; A(i + n_dash, j) = Wi;
    }else{
        //Pinned vertex, so contributions will be added to RHS i.e. B
        //Entries (i,j'), (n' + i, j'), (i, p + j'), (n' + i, p + j') will be modified
        //where j' = pinnedPointsMap[j]
        j = pinnedPointsMap.at(j);
        B(i, j) = Wr; B(i + n_dash, j + p) = Wr;
        B(i, j + p) = -Wi; B(i + n_dash, j) = Wi;
    }
}

//LSCM(Least Sqaures Conformal Mapping) Parameterization [Levy et al. 2002]
void getLSCMParameterization(const matXd &V, const vecXi &U, const matXi& F, matXd& T, int numPinnedPoints){ 
    int n_dash = F.rows();
    int n = V.rows();

    //Fix points
    //We fix on the boundary to avoid triangle flips
    vector<int> pinnedPoints = getPinnedPoints(V, U, 1);
    int p = pinnedPoints.size();

    //We design getPinnedPoints in such a way that the points are already sorted
    //in order in which they also appear in the vertex list
    //So,  i < j <=> pinnedPoints[i] < pinnedPoints[j] 

    //Map for pinned points to obtain their indices and for fast search of free points
    unordered_map<int, int> pinnedPointsMap;
    for(int i = 0; i < p; i++){
        pinnedPointsMap[pinnedPoints[i]] = i;
    }

    //Map for effective indices of free points
    unordered_map<int, int> effectiveIndexMap;
    int i, j;
    i = j = 0;
    while(i < n){
        if(pinnedPointsMap.find(i) == pinnedPointsMap.end()){
            if(i > pinnedPoints[j]){
                j++;
            }
            effectiveIndexMap[i] = i - j;   
        }
        i++;
    }

    //Solution matrices
    matXd A(2 * n_dash, 2 * (n - p)); A.setZero();
    matXd B(2 * n_dash, 2 * p); B.setZero();
    vecXd b; b.setZero(2 * n_dash);

    //Optimization for each face 
    for(unsigned int i = 0; i < n_dash; i++){
        //Obtain the coordinates of vertices in local coordinate space of triangle 
        //3D coordinates of vertices
        vec3d v0 = V.row(F(i, 0)).transpose();
        vec3d v1 = V.row(F(i, 1)).transpose();
        vec3d v2 = V.row(F(i, 2)).transpose();

        //Local coordinates of vertices
        double x0, y0, x1, y1, x2, y2;
        double s = (v1 - v0).norm();
        double t = (v2 - v0).norm();
        double cosine = (v1 - v0).dot(v2 - v0) / (s * t);

        //The actual choice of (x,y) doesn't matter as long as the lengths are angles are consistent
        x0 = y0 = 0.0f;
        x1 = s; y1 = 0.0f;
        x2 = t * cosine; y2 = t * sqrt(1 - cosine * cosine); 

        //From the area of the triangle
        double root_dT = sqrt((v1 - v0).cross(v2 - v0).norm());

        //Obtaining the weights
        //We represent the weights as (Wr^T Wi^T)^T instead of (W1 W2 W3) to split into real and imaginary parts
        //Wr = [W0r W1r W2r] and Wi = [W0i W1i W2i]
        double Wr0, Wr1, Wr2, Wi0, Wi1, Wi2;
        Wr0 = (x2 - x1) / root_dT; Wi0 = (y2 - y1) / root_dT;
        Wr1 = (x0 - x2) / root_dT; Wi1 = (y0 - y2) / root_dT;
        Wr2 = (x1 - x0) / root_dT; Wi2 = (y1 - y0) / root_dT;

        //Filling the solution matrices
        //Vertex 0
        fillMatrix(A, B, F(i,0), i, n_dash, n, p, Wr0, Wi0, pinnedPoints, pinnedPointsMap, effectiveIndexMap);
        //Vertex 1
        fillMatrix(A, B, F(i,1), i, n_dash, n, p, Wr1, Wi1, pinnedPoints, pinnedPointsMap, effectiveIndexMap);
        //Vertex 2
        fillMatrix(A, B, F(i,2), i, n_dash, n, p, Wr2, Wi2, pinnedPoints, pinnedPointsMap, effectiveIndexMap);
    }

    //Solution matrices Uf(free) and Up(pinned)
    //Note that solution vectors will be in interleaved form
    //All reals(u) followed by all imaginaries(v)

    //Solution for pinned points
    vecXd Up(2 * p);
    //Require only four points since we have two pinned points
    //We pin the points to (-1,-1) and (1,1)
    //More general algorithm pending
    Up << 0.0f, 1.0f, 0.0f, 1.0f;

    //Solution for free points
    b = -B * Up;

    //Solving the linear system
    SparseMatrix<double> A_sparse = A.sparseView();
    LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
    solver.compute(A_sparse);
    matXd Uf = solver.solve(b);

    //Accumulating the solution
    T.resize(n, 2);
    for(int i = 0; i < n; i++){
        if(pinnedPointsMap.find(i) == pinnedPointsMap.end()){
            //Free vertex
            int j = effectiveIndexMap.at(i);
            T.row(i) = vec2d(Uf(j), Uf(j + n - p)); 
        }else{
            //Pinned vertex
            int j = pinnedPointsMap.at(i);
            T.row(i) = vec2d(Up(j), Up(j + p));
        }
    }
}

double calculateEnergy(double sigma1, double sigma2, int type){
    if(type == 0){
        //Conformal Energy 
        double diff = sigma1 - sigma2;
        return 0.5 * diff * diff; 
    }else if(type == 1){
        //MIPS Energy 
        double div = sigma1 / sigma2;
        return (div + 1.0 / div);
    }
}

double getDistortion(const matXd &V, const matXd &T, const matXi &F){
    double totalDistortion = 0.0f;
    double totalArea = 0.0f; 
    //Iterate over all triangles
    for(unsigned int i = 0; i < F.rows(); i++){
        vec3d v0 = V.row(F(i, 0)).transpose();
        vec3d v1 = V.row(F(i, 1)).transpose();
        vec3d v2 = V.row(F(i, 2)).transpose();
        vec2d t0 = T.row(F(i, 0)).transpose();
        vec2d t1 = T.row(F(i, 1)).transpose();
        vec2d t2 = T.row(F(i, 2)).transpose();

        //Calculate distortion at the centroid of the triangle 
        vec3d v = (v0 + v1 + v2) / 3.0f;
        vec2d t = (t0 + t1 + t2) / 3.0f;

        //Calculate derivatives at centroid as the average of derivatives with respect to the vertices
        matXd Jt(2,3);
        //Derivative with respect to t0 = u
        vec3d dvdt0; dvdt0 << 0.0f, 0.0f, 0.0f;
        dvdt0 += (v - v0) / (t[0] - v0[0]);
        dvdt0 += (v - v1) / (t[0] - v1[0]);
        dvdt0 += (v - v2) / (t[0] - v2[0]);
        dvdt0 = dvdt0 / 3.0f;
        //Derivative with respect to t1 = v
        vec3d dvdt1; dvdt1 << 0.0f, 0.0f, 0.0f;
        dvdt1 += (v - v0) / (t[1] - v0[1]);
        dvdt1 += (v - v1) / (t[1] - v1[1]);
        dvdt1 += (v - v2) / (t[1] - v2[1]);
        dvdt1 = dvdt1 / 3.0f;
        //Assembling the Jacobian
        Jt.row(0) = dvdt0.transpose();
        Jt.row(1) = dvdt1.transpose();
        
        //The first fundamental form (size 2 x 2)
        matXd E = Jt * Jt.transpose();

        //Elements of E
        double a = E(0,0);
        double b = E(0,1);
        double c = E(1,1);

        //Calculating the singular values of the Jacobian
        double aPlusc = a + c;
        double sqrtTerm = sqrt((a-c) * (a-c) + 4 * b * b);
        double sigma1 = sqrt(0.5f * (aPlusc + sqrtTerm));
        double sigma2 = sqrt(0.5f * (aPlusc - sqrtTerm));
        
        //Obtaining the distortion using energy formulations 
        double area = ((v1 - v0).cross(v2 - v0).norm()) / 2.0f;
        totalArea += area;
        totalDistortion += area * calculateEnergy(sigma1, sigma2, 0);
    }
    return totalDistortion / totalArea;
}
