// Defines mesh class that replace the mesh toolkit in MATLAB.
// Uses DCEL data structure.

#ifndef MESH
#define MESH

#include <Eigen/Eigen>
#include <vector>

class Mesh;
class HalfEdge;
class Vertex;
class Face;

class Mesh {
    public:
        std::vector<std::shared_ptr<Vertex>> vertices;
        std::vector<std::shared_ptr<HalfEdge>> half_edges;
        std::vector<std::shared_ptr<Face>> faces;

        Mesh(std::string filename);
        
        void rotate (const int axis, const float degrees);
        int genus();
        float total_surface_area();
        float volume();

    private:

        void io_load_off(std::string filename);
        void io_load_obj(std::string filename);
        void load_vertex(Eigen::Vector3f vert);
        void load_face_with_half_edges(Eigen::Vector3i f);


        void compute_normal_per_vertex();
};

class HalfEdge{
    public:
        std::shared_ptr<HalfEdge> next; // next half edge
        std::shared_ptr<HalfEdge> twin; // twin half edge
        std::shared_ptr<Face> face; // incident face
        std::shared_ptr<Vertex> vertex; // origin vertex

        HalfEdge(std::shared_ptr<Vertex> v): vertex{v} {}
        bool has_twin(){return this->twin != nullptr;}
};

class Vertex {
    public:
        Eigen::Vector3f position;
        std::shared_ptr<HalfEdge> he; // A pointer to one of the Half Edges that starts from this Vertex 
        Eigen::Vector3f normal;
        Vertex(Eigen::Vector3f pos): position(pos){}
        std::vector<std::shared_ptr<Vertex>> neighbor_vertices();
        std::vector<std::shared_ptr<HalfEdge>> neighbor_half_edges();
};


class Face {
    public:
        std::shared_ptr<HalfEdge> he; // A pointer to one of the HalfEdges that lies on this Face 

        std::vector<std::shared_ptr<Vertex>> vertices;

        Face(std::shared_ptr<Vertex> v0, std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2){
            this->vertices.push_back(v0);
            this->vertices.push_back(v1);
            this->vertices.push_back(v2);
        }

        void compute_normal();
        std::vector<std::shared_ptr<Vertex>> get_vertices();
        float area();
        float signed_volume();
        Eigen::Vector3f normal;

};

#endif