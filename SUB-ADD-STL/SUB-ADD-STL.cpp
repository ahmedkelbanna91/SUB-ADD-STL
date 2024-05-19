#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/centroid.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/STL.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>


//using namespace CGAL;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Plane_3 Plane;
typedef std::vector<Point> Points;
typedef std::vector<std::array<std::size_t, 3>> Facets;

//// Function to compute the fitting plane for the bottom of the mesh
//Plane fit_plane_to_bottom(const Surface_mesh& mesh) {
//    std::vector<Point> points;
//    for (auto v : mesh.vertices()) {
//        points.push_back(mesh.point(v));
//    }
//    Plane plane;
//    CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), plane, CGAL::Dimension_tag<0>());
//    return plane;
//}
//
//// Function to align the mesh with the Z=0 plane
//void align_with_z_plane(Surface_mesh& mesh, const Plane& plane) {
//    Vector normal = plane.orthogonal_vector();
//    Vector up(0, 0, 1);
//
//    // Calculate rotation axis
//    Vector axis = CGAL::cross_product(normal, up);
//    double angle = std::acos(normal * up / std::sqrt(normal.squared_length()));
//
//    // Calculate rotation matrix components
//    double c = std::cos(angle);
//    double s = std::sin(angle);
//    double t = 1 - c;
//    double x = axis.x();
//    double y = axis.y();
//    double z = axis.z();
//
//    CGAL::Aff_transformation_3<Kernel> rotation(
//        t * x * x + c, t * x * y - s * z, t * x * z + s * y,
//        t * x * y + s * z, t * y * y + c, t * y * z - s * x,
//        t * x * z - s * y, t * y * z + s * x, t * z * z + c
//    );
//
//    for (auto v : mesh.vertices()) {
//        Point p = mesh.point(v);
//        mesh.point(v) = rotation(p);
//    }
//}

Point compute_centroid(const Surface_mesh& mesh) {
    std::vector<Point> vertices;
    for (auto v : mesh.vertices()) {
        vertices.push_back(mesh.point(v));
    }
    return CGAL::centroid(vertices.begin(), vertices.end());
}

void translate_mesh(Surface_mesh& mesh, const Vector& translation_vector) {
    std::cout << "      Applying translation: " << translation_vector << std::endl;
    for (auto v : mesh.vertices()) {
        mesh.point(v) = mesh.point(v) + translation_vector;
    }
}

// Utility to repair and validate a surface mesh
bool repair_and_validate_mesh(Surface_mesh mesh) {
    using namespace CGAL::Polygon_mesh_processing;

    // Remove isolated vertices
    remove_isolated_vertices(mesh);
    // fix non-manifold edges or vertices
    duplicate_non_manifold_vertices(mesh);
    // Stitch Borders:
    stitch_borders(mesh);

    // Validate the mesh after repair
    return is_valid_polygon_mesh(mesh);
}

double calculate_target_edge_length(const Surface_mesh& mesh) {
    CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    double max_dim = std::max({ bbox.xmax() - bbox.xmin(), bbox.ymax() - bbox.ymin(), bbox.zmax() - bbox.zmin() });
    return max_dim / 100.0;  // Example: Set edge length to 1% of the largest dimension
}

void remesh_surface_mesh(Surface_mesh& mesh, double target_edge_length) {
    using namespace CGAL::Polygon_mesh_processing;
    // Report the initial state of the mesh
    std::cout << "Starting remeshing..." << std::endl;
    std::cout << "Initial number of vertices: " << num_vertices(mesh) << std::endl;
    std::cout << "Initial number of edges: " << num_edges(mesh) << std::endl;
    std::cout << "Initial number of faces: " << num_faces(mesh) << std::endl;

    // Perform the remeshing
    isotropic_remeshing(
        faces(mesh),
        target_edge_length,
        mesh,
        parameters::number_of_iterations(1)  // Adjust the number of iterations as needed
    );

    // Report the final state of the mesh
    std::cout << "Remeshing completed." << std::endl;
    std::cout << "Final number of vertices: " << num_vertices(mesh) << std::endl;
    std::cout << "Final number of edges: " << num_edges(mesh) << std::endl;
    std::cout << "Final number of faces: " << num_faces(mesh) << std::endl;
}

bool read_STL(const std::string& filename, Points& points, Facets& facets) {
    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        std::cerr << "Error: Cannot open the STL file " << filename << std::endl;
        return false;
    }
    CGAL::IO::set_binary_mode(input);
    if (!CGAL::IO::read_STL(input, points, facets)) {
        std::cerr << "Error: Cannot read the STL file " << filename << std::endl;
        return false;
    }
    return true;
}

bool read_STL(const std::string& filename, Surface_mesh& mesh) {
    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        std::cerr << "Error: Cannot open the STL file to read" << filename << std::endl;
        return false;
    }
    CGAL::IO::set_binary_mode(input);
    if (!CGAL::IO::read_STL(input, mesh)) {
        std::cerr << "Error: Cannot read the STL file " << filename << std::endl;
        return false;
    }
    return true;
}

bool write_STL(const std::string& filename, const Surface_mesh& mesh) {
    std::ofstream output(filename, std::ios::binary);
    //std::cout << "Writing " << points.size() << " points and " << facets.size() << " facets to file." << std::endl;
    if (!output) {
        std::cerr << "Error: Cannot open the STL file to write" << filename << std::endl;
        return false;
    }

    if (!CGAL::is_valid_polygon_mesh(mesh)) {
        std::cerr << "Error: Mesh is not valid." << std::endl;
        // Repair and validate the mesh
        if (repair_and_validate_mesh(mesh)) {
            std::cout << "Mesh repaired." << std::endl;
        }
        else {
            std::cerr << "Error: Failed to repair or validate the mesh." << std::endl;
            return false;
        }
    }

    CGAL::IO::set_binary_mode(output);
    if (CGAL::IO::write_STL(output, mesh)) {
        output.flush();  // Flush the stream to ensure all data is written
        output.close();  // Close the file
        return true;
    }
    else {
        std::cerr << "Error: Cannot write the STL file: " << filename << std::endl;
        output.close();  // Ensure file is closed even on failure
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: SUB-ADD-STL.exe Name (fixture - tag) + model     (V1.0 CreatedByBanna)" << std::endl;
        return 1;
    }

    std::string name = argv[1];
    std::string fixture = "";
    std::string tag = "";
    std::string model = "";
    if (argc > 2) fixture = argv[2];
    if (argc > 3) tag = argv[3];
    if (argc > 4) model = argv[4];


    Surface_mesh poly1, poly2, poly3, result_subtraction, result_addition;
    Points points1, points2, points3;
    Facets facets1, facets2, facets3;
    if (!read_STL(tag, points1, facets1)) {
        std::cerr << "Error: Cannot open the Tag STL file" << std::endl;
        return EXIT_FAILURE;
    }
    if (!read_STL(fixture, points2, facets2)) {
        std::cerr << "Error: Cannot open the Fixture STL file" << std::endl;
        return EXIT_FAILURE;
    }
    if (!read_STL(model, points3, facets3)) {
        std::cerr << "Error: Cannot open the Model STL file" << std::endl;
        return EXIT_FAILURE;
    }
    // Convert STL data to Surface_mesh
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points1, facets1, poly1);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points2, facets2, poly2);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points3, facets3, poly3);


    //// Read STL files
    //if (!read_STL(tag, poly1)) {
    //    std::cerr << "Error: Cannot open the Tag STL file" << std::endl;
    //    return EXIT_FAILURE;
    //} 
    //if (!read_STL(fixture, poly2)) {
    //    std::cerr << "Error: Cannot open the Fixture STL file" << std::endl;
    //    return EXIT_FAILURE;
    //}
    //if (!read_STL(model, poly3)) {
    //    std::cerr << "Error: Cannot open the Model STL file" << std::endl;
    //    return EXIT_FAILURE;
    //}

    //double target_edge_length = calculate_target_edge_length(poly3);
    ////std::cout << "edge length = " << target_edge_length  <<std::endl;
    ////target_edge_length = 0.3;
    //remesh_surface_mesh(poly3, target_edge_length);



    //Plane plane = fit_plane_to_bottom(poly3);
    //align_with_z_plane(poly3, plane);
    // Orienting poly3 and computing its centroid
    CGAL::Polygon_mesh_processing::orient(poly3);
    Point centroid = compute_centroid(poly3);
    // Translate poly3 based on the centroid
    Kernel::Vector_3 translation(-centroid.x(), -centroid.y() + 7, 0); // Adjust Z as well if needed
    translate_mesh(poly3, translation);
   
    // Subtract poly1 from poly2
    if (!CGAL::Polygon_mesh_processing::corefine_and_compute_difference(poly2, poly1, result_subtraction)) {
        std::cerr << "Error: Subtraction operation failed." << std::endl;
        return EXIT_FAILURE;
    } 

    // Add result of subtraction to poly3
    if (!CGAL::Polygon_mesh_processing::corefine_and_compute_union(result_subtraction, poly3, result_addition)) {
        std::cerr << "Error: Addition operation failed." << std::endl;
        return EXIT_FAILURE;
    }

    // Write the result to an STL file in binary format
    if (!write_STL(name, poly3)) {
        std::cerr << "Error: Cannot write the STL file" << std::endl;
        return EXIT_FAILURE;
    } 
    else { std::cout << "      Operation completed successfully." << std::endl; }
    return EXIT_SUCCESS;
}