#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/Surface_mesh.h>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <filesystem>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef Mesh::Face_index Face_index;
typedef CGAL::Aff_transformation_3<Kernel> Transformation;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = std::filesystem;
bool DEBUG = false;


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
constexpr double deg_to_rad(double degrees) {
    return degrees * M_PI / 180.0;
}


Point compute_centroid(const Mesh& mesh, Mesh::Face_index face) {
    Point centroid(0, 0, 0);
    int num_vertices = 0;
    for (auto v : vertices_around_face(mesh.halfedge(face), mesh)) {
        centroid = CGAL::ORIGIN + (centroid - CGAL::ORIGIN) + (mesh.point(v) - CGAL::ORIGIN);
        ++num_vertices;
    }
    return CGAL::ORIGIN + (centroid - CGAL::ORIGIN) / num_vertices;
}

void flatten_face_on_z_plane(Mesh& mesh) {
    Mesh::Face_index target_face;
    double max_area = 0;
    for (auto face : mesh.faces()) {
        double area = PMP::face_area(face, mesh);
        if (area > max_area) {
            max_area = area;
            target_face = face;
        }
    }

    Vector face_normal = PMP::compute_face_normal(target_face, mesh);
    Point centroid = compute_centroid(mesh, target_face);
    Vector translation_vector = Vector(0, 0, -centroid.z());
    Transformation translation(CGAL::TRANSLATION, translation_vector);
    PMP::transform(translation, mesh);
    face_normal = PMP::compute_face_normal(target_face, mesh);

    Eigen::Vector3d axis(face_normal.x(), face_normal.y(), face_normal.z());
    Eigen::Vector3d target_axis(0, 0, -1);
    Eigen::Vector3d rotation_axis = axis.cross(target_axis).normalized();
    double angle = std::acos(axis.dot(target_axis) / (axis.norm() * target_axis.norm()));

    if (rotation_axis.norm() > 1e-6) {
        Eigen::Matrix3d rotation_matrix;
        rotation_matrix = Eigen::AngleAxisd(angle, rotation_axis);

        Transformation rotation(
            rotation_matrix(0, 0), rotation_matrix(0, 1), rotation_matrix(0, 2), 0,
            rotation_matrix(1, 0), rotation_matrix(1, 1), rotation_matrix(1, 2), 0,
            rotation_matrix(2, 0), rotation_matrix(2, 1), rotation_matrix(2, 2), 0,
            1
        );
        PMP::transform(rotation, mesh);
    }

    centroid = compute_centroid(mesh, target_face);
    translation_vector = Vector(0, 0, -centroid.z());
    Transformation final_translation(CGAL::TRANSLATION, translation_vector);
    PMP::transform(final_translation, mesh);
}

void combine_meshes(Mesh& mesh) {

    auto vpm = get(CGAL::vertex_point, mesh);
    std::vector<std::size_t> component_ids(num_faces(mesh));
    auto component_map = CGAL::make_property_map(component_ids);
    std::size_t num = PMP::connected_components(mesh, component_map, PMP::parameters::vertex_point_map(vpm));

    if (DEBUG) std::cout << "      Number of Shells: " << num << std::endl;
    if (num > 1) {
        Mesh combined, temp, temp_union;
        for (std::size_t i = 0; i < num; ++i) {
            for (Face_index f : faces(mesh)) {
                if (component_map[f] == i) {
                    auto h = halfedge(f, mesh);
                    do {
                        temp.add_vertex(mesh.point(target(h, mesh)));
                        h = next(h, mesh);
                    } while (h != halfedge(f, mesh));
                    temp.add_face();
                }
            }

            if (combined.number_of_vertices() == 0) {
                combined = temp;
            }
            else {
                PMP::corefine_and_compute_union(combined, temp, temp_union);
                combined = temp_union;
            }
        }
        mesh = combined;
    }
}

void compute_centroid(Mesh& mesh, Point& centroid) {
    if (DEBUG) std::cout << "      Computing Centroid." << std::endl;
    std::vector<Point> vertices;
    for (auto v : mesh.vertices()) {
        vertices.push_back(mesh.point(v));
    }
    centroid = CGAL::centroid(vertices.begin(), vertices.end());
}

void translate_mesh(Mesh& mesh, const Vector& translation_vector) {
    if (DEBUG) std::cout << "      Applying translation: " << translation_vector << std::endl;
    for (auto v : mesh.vertices()) {
        mesh.point(v) = mesh.point(v) + translation_vector;
    }
}

void rotate_mesh(Mesh& mesh, double x_deg, double y_deg, double z_deg) {
    double rot_x = deg_to_rad(x_deg);
    double rot_y = deg_to_rad(y_deg);
    double rot_z = deg_to_rad(z_deg);

    double cos_x = std::cos(rot_x), sin_x = std::sin(rot_x);
    Transformation rot_mtx_x(
        1, 0, 0, 0,
        0, cos_x, -sin_x, 0,
        0, sin_x, cos_x, 0,
        1
    );

    double cos_y = std::cos(rot_y), sin_y = std::sin(rot_y);
    Transformation rot_mtx_y(
        cos_y, 0, sin_y, 0,
        0, 1, 0, 0,
        -sin_y, 0, cos_y, 0,
        1
    );

    double cos_z = std::cos(rot_z), sin_z = std::sin(rot_z);
    Transformation rot_mtx_z(
        cos_z, -sin_z, 0, 0,
        sin_z, cos_z, 0, 0,
        0, 0, 1, 0,
        1
    );

    Transformation combined = rot_mtx_x * rot_mtx_y * rot_mtx_z;
    PMP::transform(combined, mesh);
}

bool repair_and_validate_mesh(Mesh& mesh) {
    PMP::remove_isolated_vertices(mesh);
    PMP::duplicate_non_manifold_vertices(mesh);
    PMP::stitch_borders(mesh);
    return CGAL::is_valid_polygon_mesh(mesh);
}

bool read_STL(const std::string& filename, Mesh& mesh) {
    fs::path filepath(filename);
    if (DEBUG) std::cout << "      Reading STL file: " << filepath.filename() << std::endl;
    if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
        std::cerr << "Error: Cannot read the STL file " << filepath.filename() << std::endl;
        return false;
    }
    return true;
}

bool write_STL(const std::string& filename, const Mesh& mesh) {
    fs::path filepath(filename);
    if (DEBUG) std::cout << "      Writting STL file." << filepath.filename() << std::endl;
    if (!CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(17))) {
        std::cerr << "Error: Cannot write the STL file: " << filepath.filename() << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-D") {
            DEBUG = true; 
            continue;
        }
        if (i + 1 < argc) {
            args[argv[i]] = argv[i + 1]; 
            i++;
        }
        else {
            std::cerr << "Missing value for " << argv[i] << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (args.find("-F") == args.end() || args.find("-T") == args.end() || args.find("-O") == args.end()) {
        std::cerr << "Usage: SUB-ADD-STL.exe -F fixture.stl -T tag.stl [-I model.stl] -O out.stl [-D]" << std::endl;
        return EXIT_FAILURE;
    }

    std::string fixture = args["-F"], tag = args["-T"], model = args["-I"], output = args["-O"];
    Mesh tag_mesh, fixture_mesh, model_mesh, subtraction_result, addition_result, result;
    Point centroid;

    if (!read_STL(tag, tag_mesh) || !read_STL(fixture, fixture_mesh)) {
        return EXIT_FAILURE;
    }

    if (!PMP::corefine_and_compute_difference(fixture_mesh, tag_mesh, subtraction_result)) {
        std::cerr << "Error: Subtraction operation failed." << std::endl;
        return EXIT_FAILURE;
    }
    //result = subtraction_result;

    if (!model.empty()) {
        if (!read_STL(model, model_mesh)) {
            return EXIT_FAILURE;
        }

        flatten_face_on_z_plane(model_mesh);
        result = model_mesh;


        //rotate_mesh (model_mesh, 0, 0, 180);
        /*settle_mesh (model_mesh);
        result = model_mesh;*/

        //PMP::orient(model_mesh);
        //compute_centroid(model_mesh, centroid);
        //translate_mesh(model_mesh, Kernel::Vector_3(-centroid.x(), -centroid.y() + 7, 0));

        //CGAL::copy_face_graph(model_mesh, addition_result);
        //CGAL::copy_face_graph(subtraction_result, addition_result);

        //combine_meshes(addition_result);
        //if (!PMP::corefine_and_compute_union(model_mesh, subtraction_result, addition_result)) {
        //    std::cerr << "Error: Addition operation failed." << std::endl;
        //    return EXIT_FAILURE;
        //}

        //result = addition_result;
    }

    if (!CGAL::is_valid_polygon_mesh(result)) {
        std::cerr << "Error: Mesh is not valid." << std::endl;
        if (repair_and_validate_mesh(result)) {
            if (DEBUG) std::cout << "Mesh repaired." << std::endl;
        }
        else {
            std::cerr << "Error: Failed to repair or validate the mesh." << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    if (!write_STL(output, result)) {
        return EXIT_FAILURE;
    }

    std::cout << "      Operation completed successfully." << std::endl;
    return EXIT_SUCCESS;
}