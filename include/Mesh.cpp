// Defines mesh helper functions that replace the mesh toolkit in MATLAB.
// Uses DCEL data structure.

#include "Mesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

// Mesh::load_vertex() - Takes position of the vertex and push it into the list
void Mesh::load_vertex(Eigen::Vector3f vert){
    std::shared_ptr<Vertex> v_ptr (new Vertex(vert));
    v_ptr->he = nullptr;
    this->vertices.push_back(v_ptr);
}

// Mesh::load_face_with_half_edges() - Take three vertices of the triangle 
// and push the relavent face and half edges into the list
void Mesh::load_face_with_half_edges(Eigen::Vector3i f) {
    std::shared_ptr<HalfEdge> e1 (new HalfEdge(vertices.at(f.x())));
    std::shared_ptr<HalfEdge> e2 (new HalfEdge(vertices.at(f.y())));
    std::shared_ptr<HalfEdge> e3 (new HalfEdge(vertices.at(f.z())));

    std::shared_ptr<Face> fa (new Face(vertices.at(f.x()), vertices.at(f.y()),vertices.at(f.z())));

    if (vertices.at(f.x())->he == nullptr) vertices.at(f.x())->he = e1;
    if (vertices.at(f.y())->he == nullptr) vertices.at(f.y())->he = e2;
    if (vertices.at(f.z())->he == nullptr) vertices.at(f.z())->he = e3;
           
    e1->next = e2;
    e2->next = e3;
    e3->next = e1;
    e1->vertex = vertices.at(f.x());
    e2->vertex = vertices.at(f.y());
    e3->vertex = vertices.at(f.z());
    e1->face = fa;
    e2->face = fa;
    e3->face = fa;

    for (auto& t : this->half_edges){
        if (t->vertex == e1->next->vertex){
            if (e1->vertex == t->next->vertex){
                e1->twin = t;
                t->twin = e1;
            }
        }
        if (t->vertex == e2->next->vertex){
            if (e2->vertex == t->next->vertex){
                e2->twin = t;
                t->twin = e2;
            }
                }
        if (t->vertex == e3->next->vertex){
            if (e3->vertex == t->next->vertex){
                e3->twin = t;
                t->twin = e3;
            }
        }
    }

    fa->he = e1;
    half_edges.push_back(e1);
    half_edges.push_back(e2);
    half_edges.push_back(e3);
    faces.push_back(fa);
}


// Iterate through all neighbour vertices around the vertex
std::vector<std::shared_ptr<Vertex>> Vertex::neighbor_vertices() {
    std::vector<std::shared_ptr<Vertex>> neighborhood;
    auto he = this->he;
    do {
        neighborhood.push_back(he->twin->vertex);
        he = he->twin->next;
    }
    while(he != this->he);
    return neighborhood; 
}


// Iterate through all half edges pointing away from the vertex
std::vector<std::shared_ptr<HalfEdge>> Vertex::neighbor_half_edges() {
    std::vector<std::shared_ptr<HalfEdge>> neighborhood;
    auto he = this->he;
    do {
        neighborhood.push_back(he);
        he = he->twin->next;
    }
    while(he != this->he);
    return neighborhood;
}

// Iterate through all member vertices of the face
std::vector<std::shared_ptr<Vertex>> Face::get_vertices() {
    std::vector<std::shared_ptr<Vertex>> member_vertices;
    auto he = this->he;
    do {
        member_vertices.push_back(he->vertex);
        he = he->next;
    }
    while(he != this->he);
    return member_vertices;
}


// Compute the area of the triangular face 
float Face::area(){
    auto vertices = Face::get_vertices();
    Eigen::Vector3f v0 = vertices.at(0)->position, v1 = vertices.at(1)->position, v2 = vertices.at(2)->position;

    float a = sqrt(pow(v0.x() - v1.x(), 2) + pow(v0.y() - v1.y(), 2) + pow(v0.z() - v1.z(), 2));
    float b = sqrt(pow(v2.x() - v1.x(), 2) + pow(v2.y() - v1.y(), 2) + pow(v2.z() - v1.z(), 2));
    float c = sqrt(pow(v0.x() - v2.x(), 2) + pow(v0.y() - v2.y(), 2) + pow(v0.z() - v2.z(), 2));
    float s = (a + b + c) / 2.f;
    
    return sqrt(s * (s - a) * (s - b) * (s - c));
}

// Compute the signed volume of the triangular face 
// reference: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf eq.(5)
float Face::signed_volume(){
    auto vertices = Face::get_vertices();
    Eigen::Vector3f v0 = vertices.at(0)->position, v1 = vertices.at(1)->position, v2 = vertices.at(2)->position;
    return (-v2.x() * v1.y() * v0.z() + v1.x() * v2.y() * v0.z() + v2.x() * v0.y() * v1.z()
            -v0.x() * v2.y() * v1.z() - v1.x() * v0.y() * v2.z() + v0.x() * v1.y() * v2.z()) / 6.f; 

}

// Use the half-edge data structure, compute the surface normal vector for each triangle face. 
void Face::compute_normal(){
    auto vertices = Face::get_vertices();
    Eigen::Vector3f v0 = vertices.at(0)->position, v1 = vertices.at(1)->position, v2 = vertices.at(2)->position;
    normal = (v2 - v1).cross(v0 - v1).normalized();
}

// Compute the genus number 
int Mesh::genus() {
    int genus = 0;
    genus = 1 - (vertices.size() - half_edges.size() / 2 + faces.size()) / 2;
    return genus;
}


// Compute the surface area of the mesh
float Mesh::total_surface_area() {
    float total_surface_area = 0;
    for (auto& f: this->faces){
        total_surface_area += f->area();
    }
    return total_surface_area;
}


// Compute the volume of the mesh 
float Mesh::volume() {
    float total_volume = 0;
    for (auto& f: this->faces){
        total_volume += f->signed_volume();
    }
    return std::fabs(total_volume);
}

// Compute the normal vector for vertex. 
void Mesh::compute_normal_per_vertex() {


    for (auto& f : this->faces) f->compute_normal();
    for (auto& v : this->vertices) {
        if (v->he == nullptr) {
            v->normal = Eigen::Vector3f::Zero();
            continue;
        }
        v->normal = Eigen::Vector3f::Zero();
        auto neighbors = v->neighbor_half_edges();
        std::vector<float> angles;
        float total = 0;
        for (int i = 0; i < neighbors.size() - 1; i++){
            angles.push_back(acos((neighbors.at(i)->vertex->position - v->position).normalized()
                    .dot((neighbors.at(i + 1)->vertex->position - v->position).normalized())));
        }
        angles.push_back(acos((neighbors.at(neighbors.size() - 1)->vertex->position - v->position).normalized()
                    .dot((neighbors.at(0)->vertex->position - v->position).normalized())));
        for (int i = 0; i < neighbors.size(); i++){
            v->normal += neighbors.at(i)->face->normal * angles.at(i);
        }
        v->normal = v->normal.normalized();
    }
}

void Mesh::io_load_off(std::string filename){
//  M = LOAD_OFF(filename) Loads a mesh from the given OFF file.
//
//  The file is expected to begin directly after the header, without
//  containing black or commented lines.
//  The function loads a triangle mesh, so the file is expected to contain
//  only triangles.

    
    // Try to open the file and check if it has benn successfully opened
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
    }
    else {
        // Check the file header to ensure the given file is an OFF
        std::string line;
        std::getline(file, line);
        if (line != "OFF") 
            std::cerr << "File " + filename + " is not a valid OFF file." << "\n";
        else{
            int v_count, f_count;

            while (std::getline(file, line)) {
                std::istringstream iss(line);
                if (line == "") continue;
                std::string prefix;
                iss >> prefix;
                if ("#" == prefix) continue;
                else {
                    v_count = std::stoi(prefix);
                    iss >> f_count;
                    break;
                }
            }

            while (v_count > 0) {
                std::getline(file, line);
                if (line == "") continue;
                std::istringstream iss(line);
                std::string prefix;
                iss >> prefix;
                float x, y, z;

                if ("#" == prefix) continue;
                else {
                    x = std::stof(prefix);
                    iss >> y >> z;
                    load_vertex(Eigen::Vector3f(x, y, z));
                    v_count--;
                }
            }

            std::cout << vertices.size() << " vertices loaded.\n";

            while (f_count > 0) {
                std::getline(file, line);
                if (line == "") continue;
                std::istringstream iss(line);
                std::string prefix;
                iss >> prefix;
                if ("#" == prefix) continue;
                else {
                    const int num_of_vertices = std::stoi(prefix);
                    std::vector<int> vs;
                    int v;
                    for (int i = 0; i < num_of_vertices; i++){
                        iss >> v;
                        vs.push_back(v);
                    }

                    for (int i = 0; i < num_of_vertices - 2; i++){
                        load_face_with_half_edges(Eigen::Vector3i(vs.at(0), vs.at(i + 1), vs.at(i + 2)));
                    }
                    f_count--;
                }
            }
            std::cout << faces.size() << " faces and " << half_edges.size() << " half-edges loaded.\n";
        }
    }
}

void Mesh::io_load_obj(std::string filename){
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
    }

    else {
        std::cout << "Now loading mesh...\n";
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string prefix;
            std::vector<int> nodes;
            iss >> prefix;

            if ("v" == prefix) {
                Eigen::Vector3f vert;
                iss >> vert[0] >> vert[1] >> vert[2];
                load_vertex(vert);
            } 

            else if ("f" == prefix) {
                std::string sub;
                int i = 3;
                while (i--){
                    iss >> sub;
                    std::string num = "";
                    for (char& c: sub) {
                        if (c == '/') break;
                        num = num + c;
                    }
                    nodes.push_back(std::stoi(num));
                }
                sub = "";
                Eigen::Vector3i face(nodes.at(0), nodes.at(1), nodes.at(2));
                face = face.array() - 1;
                load_face_with_half_edges(face);

                iss >> sub;
                if (sub != "") { // Assume it's a quad mesh in this case, divide it into two triangles
                    std::string num = "";
                    for (char& c: sub) {
                        if (c == '/') break;
                        num = num + c;
                    }
                    nodes.push_back(std::stoi(num));
                    Eigen::Vector3i face1(nodes.at(0), nodes.at(2), nodes.at(3));
                    face1 = face1.array() - 1;
                    load_face_with_half_edges(face1);
                }
            }
        }
    }
}

Mesh::Mesh(std::string filename){
//  M = INIT(filename) Loads a triangle mesh from the given file.
//
    const std::string file_extension = filename.substr(filename.size() - 4, 4);
    if (file_extension == ".off") io_load_off(filename);
}