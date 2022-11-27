#pragma once

#include "common.hpp"

struct Vertex;
struct Edge; 
struct Face;

struct Vertex {
    int id;
    vec3d pos;
    Edge* edge;
};

struct Edge {
    int id;
    Vertex* startVertex;
    Face* face;
    Edge* next;
    Edge* prev;
    Edge* twin;
    Edge(int id, Vertex* startVertex, Face* face){
        this->id = id;
        this->startVertex = startVertex;
        this->face = face;
        this->next = NULL;
        this->prev = NULL;
        this->twin = NULL;
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

void Graph::getBoundary(vecXi& Vb){
    //We assume that this is invoked only on meshes with a boundary i.e. boundary exists
    //Find a half edge with no face mapped no it
    Vertex* firstVertex;
    for(unsigned int i = 0; i < edges.size(); i++){
        if(edges[i]->face == NULL){
            firstVertex = edges[i]->startVertex;
            break;
        }
    }
    Vertex* currentVertex = firstVertex;
    vector<int> boundary;
    while(true){
        //Insert the current vertex into the boundary
        boundary.push_back(currentVertex->id);
        //Find the next vertex
        Edge* currentEdge = currentVertex->edge;
        while(currentEdge->twin->face){
            currentEdge = currentEdge->twin->next;
        }
        currentVertex = currentEdge->twin->startVertex;
        if(currentVertex == firstVertex){
            break;
        }
    }
    Vb.resize(boundary.size());
    for(unsigned int i = 0; i < boundary.size(); i++){
        Vb(i) = boundary[i];
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