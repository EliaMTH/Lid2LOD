#include <cinolib/meshes/meshes.h>
// #include <cinolib/gl/glcanvas.h>
// #include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/profiler.h>

#include "../src/ground.h"
#include "../src/buildings.h"
#include "../src/city.h"
#include "../src/auxiliary.h"

using namespace cinolib;

int main(int argc, char *argv[])
{
    if (argc != 6) {
        std::cout << "Welcome to TriCity! Please specify the followings:\n"
                     "- the path to the boundary .off file\n"
                     "- the path to the streets .csv file\n"
                     "- the path to the folder containing the triangulated buildings .off files\n"
                     "- the path to the output folder\n"
                     "- if you want to visualize the output (1) or not (0)\n";
        exit(0);
    }
    std::string boundary_path  = argv[1];
    std::string streets_path   = argv[2];
    std::string buildings_path = argv[3];
    std::string output_path    = argv[4];
    bool VISUALIZE             = atoi(argv[5]);

    /*************** CREATE GROUND AND CITY MESH ****************/

    Profiler prof;
    prof.push("Total");

    prof.push("Load data");
    Polygonmesh<> boundary(boundary_path.c_str());
    std::vector<std::string> buildings_dirs = load_buildings_dirs(buildings_path);
    std::vector<std::tuple<int,vec3d,vec3d>> streets = load_streets(streets_path);
    std::string msg = " found " + std::to_string(buildings_dirs.size()) + " buildings and "
                      + std::to_string(streets.size()) + " streets.";
    prof.pop(true, msg);

    bool SAVE = !output_path.empty();
    if (SAVE) open_directory(output_path);

    prof.push("Create ground mesh");
    Trimesh<> ground = create_ground_mesh(boundary, buildings_dirs, streets);
    prof.pop();
    if (SAVE) {
        ground.save((output_path + "/ground.obj").c_str());
    }

    prof.push("Create buildings mesh");
    Trimesh<> buildings = create_buildings_mesh(buildings_dirs);
    prof.pop();
    if (SAVE) {
        buildings.save((output_path + "/buildings.obj").c_str());
    }

    prof.push("Create city mesh");
    Trimesh<> city = create_city_mesh(ground, buildings);
    prof.pop();
    if (SAVE) {
        city.save((output_path + "/city.obj").c_str());
        // print_poly_labels(city, output_path + "/polys_buildings_IDs.csv");
        // print_edge_labels(city, output_path + "/edges_streets_IDs.csv");
    }

    prof.pop();

    /*************** DISPLAY THE RESULT ****************/

    // if (VISUALIZE) {
    //     DrawablePolygonmesh<> DM = convert_to_drawable(ground);
    //     DM.translate(-DM.centroid());
    //     DM.update_bbox();
    //     // DM.edge_mark_boundaries();
    //     DM.poly_color_wrt_label();
    //     for (uint pid=0; pid<DM.num_polys(); ++pid) {
    //         if (DM.poly_data(pid).label == -1) {
    //             DM.poly_data(pid).color = Color::WHITE();
    //         }
    //     }
    //     DM.updateGL();

    //     GLcanvas gui(1000, 1000);
    //     gui.push(&DM);

    //     // #include "../src/streets.h"
    //     // std::vector<std::tuple<int,vec3d,vec3d>> streets;
    //     // load_streets(streets_path, streets);
    //     // DrawableSegmentSoup streets_soup = visualize_streets(streets, ground_boundary.centroid());
    //     // gui.push(&streets_soup);

    //     SurfaceMeshControls<DrawablePolygonmesh<>> menu(&DM, &gui, "City");
    //     gui.push(&menu);
    //     return gui.launch();
    // }
    return 0;
}
