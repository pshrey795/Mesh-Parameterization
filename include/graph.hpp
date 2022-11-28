#pragma once

#include "parameterization.hpp"

struct Vertex;
struct Edge; 
struct Face;

struct Vertex {
    int id;
    vec3d pos;
    Edge* edge;
    bool isBoundary; 
};

struct Edge {
    int id;
    Vertex* startVertex;
    Face* face;
    Edge* next;
    Edge* prev;
    Edge* twin;
    bool isBoundary;
    Edge(int id, Vertex* startVertex, Face* face){
        this->id = id;
        this->startVertex = startVertex;
        this->face = face;
        this->next = NULL;
        this->prev = NULL;
        this->twin = NULL;
        this->isBoundary = false;
    }
};

struct Face {
    int id;
    Edge* edge; 
    int indices[3];
    Face(int id, int a, int b, int c){
        this->id = id;
        indices[0] = a;
        indices[1] = b;
        indices[2] = c;
        this->edge = NULL;
    }
};

class Graph {
    public:
        vector<Vertex*> vertices;
        vector<Edge*> edges;
        vector<Face*> faces;
        vector<vector<int>> adjList;

        Graph(const matXd& V, const matXi& F);
        void getBoundary(vecXi &Vb);
        vector<unsigned int> boundaryVertices;

        void getInternalParameterization(const matXd& Tb, matXd& T, int method);
        vector<int> getNeighbours(int i);

    private:
        double getWeight(int i, int j, int method);

};

Graph::Graph(const matXd& V, const matXi& F){
    unsigned int n = V.rows();
    unsigned int m = F.rows();

    //Preparing adjacency list
    for(unsigned int i = 0; i < n; i++){
        vector<int> adj;
        for(unsigned int j = 0; j < n; j++){
            adj.push_back(-1);
        }
        adjList.push_back(adj);
    }

    //Adding vertices
    for(unsigned int i = 0; i < n; i++){
        Vertex* v = new Vertex;
        v->id = i;
        v->pos = V.row(i).transpose();
        v->edge = NULL;
        v->isBoundary = false;
        vertices.push_back(v);
    }

    unsigned int i = 0;
    while(i < m){
        int a, b, c;

        //Vertex indices
        a = F(i, 0);
        b = F(i, 1);
        c = F(i, 2);

        //Create face and half edges
        Face* face = new Face(i, a, b, c);
        Edge* e1 = new Edge(3*i, vertices[a], face);
        Edge* e2 = new Edge(3*i+1, vertices[b], face);
        Edge* e3 = new Edge(3*i+2, vertices[c], face);

        //Linking edges to vertices if not done already
        if(vertices[a]->edge == NULL){
            vertices[a]->edge = e1;
        }
        if(vertices[b]->edge == NULL){
            vertices[b]->edge = e2;
        }
        if(vertices[c]->edge == NULL){
            vertices[c]->edge = e3;
        }

        //Log twin edges
        adjList[a][b] = 3*i;
        adjList[b][c] = 3*i+1;
        adjList[c][a] = 3*i+2;

        //Is the order of the edges important? (Clockwise or Anti-Clockwise)
        e1->next = e2;
        e2->next = e3;
        e3->next = e1;
        e1->prev = e3;
        e2->prev = e1;
        e3->prev = e2;

        //Linking edges to face
        face->edge = e1;

        //Adding new face and edges to the graph
        faces.push_back(face);
        edges.push_back(e1);
        edges.push_back(e2);
        edges.push_back(e3);

        //Next iteration
        i += 1;
    }
    int currentSize = edges.size();

    //Linking twin edges
    for(unsigned int a = 0; a < n; a++){
        for(unsigned int b = 0; b < n; b++){
            if(adjList[a][b] != -1){
                if(adjList[b][a] != -1){
                    //The other way around will be handled in its own iteration 
                    edges[adjList[a][b]]->twin = edges[adjList[b][a]]; 
                    edges[adjList[b][a]]->twin = edges[adjList[a][b]];
                }
            }else{
                if(adjList[b][a] != -1){
                    Vertex* newStartVertex = edges[adjList[b][a]]->next->startVertex;
                    adjList[a][b] = currentSize++;
                    Edge* newEdge = new Edge(adjList[a][b], newStartVertex, NULL);
                    newEdge->twin = edges[adjList[b][a]];
                    edges[adjList[b][a]]->twin = newEdge;
                    edges.push_back(newEdge);
                }
            }
        }
    }
}

//Debug: done
void Graph::getBoundary(vecXi& Vb){
    //We assume that this is invoked only on meshes with a boundary i.e. boundary exists
    //Find a half edge with no face mapped no it
    Vertex* firstVertex;
    for(unsigned int i = 0; i < edges.size(); i++){
        if(edges[i]->face == NULL){
            edges[i]->isBoundary = true;
            firstVertex = edges[i]->startVertex;
            break;
        }
    }
    Vertex* currentVertex = firstVertex;
    while(true){
        //Insert the current vertex into the boundary
        currentVertex->isBoundary = true;
        boundaryVertices.push_back(currentVertex->id);
        //Find the next vertex
        Edge* currentEdge = currentVertex->edge;
        while(currentEdge->twin->face){
            currentEdge = currentEdge->twin->next;
        }
        currentVertex = currentEdge->twin->startVertex;
        if(currentVertex == firstVertex){
            break;
        }else{
            currentEdge->twin->isBoundary = true;
        }
    }
    Vb.resize(boundaryVertices.size());
    for(unsigned int i = 0; i < boundaryVertices.size(); i++){
        Vb(i) = boundaryVertices[i];
    }
}

void getMeshFromGraph(Graph& graph, matXd& V, matXi& F){

    unsigned int n = graph.vertices.size();
    unsigned int m = graph.faces.size();

    V.resize(n, 3);
    F.resize(m, 3);

    for(unsigned int i = 0; i < n; i++){
        V.row(i) = graph.vertices[i]->pos.transpose();
    }

    for(unsigned int i = 0; i < m; i++){
        F(i, 0) = graph.faces[i]->indices[0];
        F(i, 1) = graph.faces[i]->indices[1];
        F(i, 2) = graph.faces[i]->indices[2];
    }
}

//For obtaining the parameter points of internal mesh vertices using affine combinations 
//Method determines the weights of the affine combination 
//Method = 0: Uniform weights(Tutte's embedding)
//Method = 1: Harmonic/Cotangent weights
//Method = 2: Mean value weights(Floater) 
void Graph::getInternalParameterization(const matXd& Tb, matXd& T, int method){
    //Given: Boundary parameterization Tb
    //Find: Complete parameterization T
    
    //Debug: Done
    //First identify the indices of the internal vertices
    vector<int> internalVertices;
    for(unsigned int i = 0; i < vertices.size(); i++){
        if(!vertices[i]->isBoundary){
            internalVertices.push_back(i);
        }
    }

    //Debug: Done
    //Map for obtaining the boundary coordinates using their indices
    unordered_map<int, vec2d> boundaryMap;
    for(unsigned int i = 0; i < boundaryVertices.size(); i++){
        boundaryMap[boundaryVertices[i]] = Tb.row(i).transpose();
    }

    //Debug: Done
    //Map for obtaining the index of the internal coordinates for querying the solution matrices
    unordered_map<int, int> internalMap;
    for(unsigned int i = 0; i < internalVertices.size(); i++){
        internalMap[internalVertices[i]] = i;
    }

    //Debug: Done
    //Preparing the linear system Au = b_u and Av = b_v
    unsigned int n = internalVertices.size();
    //LHS
    matXd A(n, n);
    A.setZero(); 
    //RHS
    vecXd b_u(n);
    vecXd b_v(n);
    b_u.setZero(); b_v.setZero();
    //Solution
    vecXd u(n);
    vecXd v(n);
    u.setZero(); v.setZero();
    for(unsigned int i = 0; i < n; i++){
        int currentIdx = internalVertices[i];
        vector<int> neighbors = this->getNeighbours(currentIdx);
        int degree = neighbors.size();
        A(i, i) = 1.0; //Diagonal entries are 1
        for(unsigned int j = 0; j < degree; j++){
            int neighborIdx = neighbors[j];
            double weight = this->getWeight(currentIdx, neighborIdx, method);
            if(true){
                weight = weight/(double)degree;
            }
            if(this->vertices[neighborIdx]->isBoundary){
                //Boundary vertex, so add contribution to RHS
                vec2d rhsVal = weight * boundaryMap[neighborIdx];
                b_u(i) += rhsVal[0];
                b_v(i) += rhsVal[1];
            }else{
                //Internal vertex, so add contribution to LHS
                neighborIdx = internalMap[neighborIdx];
                A(i, neighborIdx) = (-1.0) * weight;                
            }
        }
    }

    //Debug: Done
    //Solving the linear system
    SparseMatrix<double> A_sparse = A.sparseView();
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(A_sparse);
    if(solver.info() != Eigen::Success){
        cout << "Eigen factorization failed!" << endl; 
    }
    u = solver.solve(b_u);
    if(solver.info() != Eigen::Success){
        cout << "Eigen solve failed for u!" << endl; 
    }
    v = solver.solve(b_v);
    if(solver.info() != Eigen::Success){
        cout << "Eigen solve failed for v!" << endl; 
    }

    //Debug: Done
    //Accumulating the solution
    T.resize(vertices.size(), 2);
    for(unsigned int i = 0; i < vertices.size(); i++){
        if(vertices[i]->isBoundary){
            T.row(i) = boundaryMap[i];
        }else{
            int idx = internalMap[i];
            T(i, 0) = u(idx);
            T(i, 1) = v(idx);
        }
    }
}

double Graph::getWeight(int i, int j, int method){
    if(method == 0){
        //Uniform weights
        return 1.0;
    }else if(method == 1){
        return 1.0;
    }else{
        return 1.0;
    }
}

//Debug: done
vector<int> Graph::getNeighbours(int i){
    vector<int> res;
    Edge* currentEdge = this->vertices[i]->edge;
    res.push_back(currentEdge->twin->startVertex->id);
    if(currentEdge->prev != NULL){
        Edge* nextRightEdge = currentEdge->prev->twin;
        while(nextRightEdge != currentEdge){
            res.push_back(nextRightEdge->twin->startVertex->id);
            if(nextRightEdge->prev != NULL){
                nextRightEdge = nextRightEdge->prev->twin;
            }else{
                break;
            }
        }
        if(nextRightEdge != currentEdge){
            Edge* nextLeftEdge = currentEdge->twin->next;
            while(nextLeftEdge != NULL){
                res.push_back(nextLeftEdge->twin->startVertex->id);
                nextLeftEdge = nextLeftEdge->twin->next;
            }
        }
    }else{
        Edge* nextLeftEdge = currentEdge->twin->next;
        if(nextLeftEdge != NULL){
            while(nextLeftEdge != currentEdge){
                if(nextLeftEdge != NULL){
                    res.push_back(nextLeftEdge->twin->startVertex->id);
                    nextLeftEdge = nextLeftEdge->twin->next;
                }else{
                    break;
                }
            }
        }
    }
    return res; 
}