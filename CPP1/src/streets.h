#ifndef STREETS
#define STREETS

// #include <cinolib/drawable_segment_soup.h>

#include "auxiliary.h"

using namespace cinolib;

std::vector<std::tuple<int,vec3d,vec3d>> load_streets(const std::string &path)
{
    std::vector<std::tuple<int,vec3d,vec3d>> streets;

    // load the streets graph
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Failed to open streets_file.\n";
        return streets;
    }
    // Skip header line
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        uint eid;
        double x0, y0, z0, x1, y1, z1;
        std::getline(ss, token, ','); eid = std::stoul(token);
        std::getline(ss, token, ','); x0 = std::stod(token);
        std::getline(ss, token, ','); y0 = std::stod(token);
        std::getline(ss, token, ','); z0 = std::stod(token);
        std::getline(ss, token, ','); x1 = std::stod(token);
        std::getline(ss, token, ','); y1 = std::stod(token);
        std::getline(ss, token, ','); z1 = std::stod(token);

        vec3d v0(x0, y0, z0);
        vec3d v1(x1, y1, z1);

        // Check if v0 = v1
        // if (v0 == v1) {
        //     std::cout << "load_streets - ERROR: the street points are the same: "
        //               << "v0 = " << v0 << ", v1 = " << v1 << "; "
        //               << "eid = " << eid << std::endl;
        //     continue;
        // }

        // // Check if the street is already in the list
        // bool found = false;
        // for (const auto &street : streets) {
        //     if ((std::get<1>(street) == v0 && std::get<2>(street) == v1) ||
        //         (std::get<1>(street) == v1 && std::get<2>(street) == v0)) {
        //         found = true;
        //         break;
        //     }
        // }
        // if (found) {
        //     std::cout << "load_streets - ERROR: the street is already in the list: "
        //               << "v0 = " << v0 << ", v1 = " << v1 << "; "
        //               << "eid = " << eid << std::endl;
        //     continue;
        // }

        // Add the street to the list
        streets.emplace_back(std::make_tuple(eid, v0, v1));
    }
    return streets;
}

// DrawableSegmentSoup visualize_streets(const std::vector<std::tuple<int,vec3d,vec3d>> &streets, const vec3d &center)
// {
//     DrawableSegmentSoup soup;
//     for (const auto &street : streets) {
//         int sid  = std::get<0>(street);
//         vec3d v0 = std::get<1>(street) - center;
//         vec3d v1 = std::get<2>(street) - center;
//         soup.push_seg(v0, v1, Color::BLUE());
//     }

//     soup.default_color = Color::BLUE();
//     soup.thickness = 1.0;
//     soup.no_depth_test = true;
//     soup.use_gl_lines = true;

//     return soup;
// }

void translate_and_project_street(std::tuple<int,vec3d,vec3d> &street,
                                  const vec3d &translation,
                                  std::unordered_map<uint, double> &z_coords)
{
    // translate the street to the desired center and project onto the z=0 plane
    std::get<1>(street) += -translation;
    std::get<2>(street) += -translation;

    uint offset = z_coords.size();
    uint id0 = offset;
    uint id1 = offset+1;
    project(id0, std::get<1>(street), z_coords);
    project(id1, std::get<2>(street), z_coords);
}

template<class M, class V, class E, class P>
double get_avg_z(const AbstractPolygonMesh<M,V,E,P> &m,
                 const uint &vid)
{
    double avg_z = 0.;
    for (uint nbr : m.adj_v2v(vid)) {
        avg_z += m.vert(nbr).z();
    }
    avg_z /= m.adj_v2v(vid).size();
    return avg_z;
}

bool is_street_inside(const std::tuple<int,vec3d,vec3d> &street,
                      const Polygonmesh<> &ground,
                      const uint pid0)
{
    vec3d v0 = std::get<1>(street);
    vec3d v1 = std::get<2>(street);

    // check that the street is inside the boundary
    if (!point_in_polygon(ground, v0, pid0) || !point_in_polygon(ground, v1, pid0)) {
        return false;
    }
    // check that the street is not crossing a building
    for (uint pid=1; pid<ground.num_polys(); ++pid) {
        if (point_in_polygon(ground, v0, pid) || point_in_polygon(ground, v1, pid)) {
            return false;
        }
    }
    return true;
}

void add_street_to_mesh(const std::tuple<int,vec3d,vec3d> &street,
                        Polygonmesh<> &m)
{
    // add the street points and edges to the mesh
    int sid  = std::get<0>(street);
    vec3d v0 = std::get<1>(street);
    vec3d v1 = std::get<2>(street);

    auto it0 = std::find(m.vector_verts().begin(), m.vector_verts().end(), v0);
    uint vid0 = it0 != m.vector_verts().end() ?
                    std::distance(m.vector_verts().begin(), it0) : m.vert_add(v0);

    auto it1 = std::find(m.vector_verts().begin(), m.vector_verts().end(), v1);
    uint vid1 = it1 != m.vector_verts().end() ?
                    std::distance(m.vector_verts().begin(), it1) : m.vert_add(v1);

    if (vid0 == vid1) {
        std::cout << "add_street_to_mesh - ERROR: the street extremities coincide: "
                  << "vid0 = " << vid0 << ", vid1 = " << vid1 << "; "
                  << v0 << " and " << v1 << std::endl;
        // assert(false);
        return;
    }

    int eid = m.edge_id(vid0, vid1);
    if (eid == -1) {
        eid = m.edge_add(vid0, vid1);
    }
    m.edge_data(eid).label = sid;
}

template<class M, class V, class E, class P>
void mark_streets_edges_and_verts(const AbstractPolygonMesh<M,V,E,P> &m_old,
                                  AbstractPolygonMesh<M,V,E,P> &m_new)
{
    // mark the edges corresponding to streets
    for (uint eid=0; eid<m_old.num_edges(); ++eid) {
        auto vids = m_old.edge_vert_ids(eid);
        int eid_tri = m_new.edge_id(vids);
        if (eid_tri == -1) {
            // std::cout << "WARNING: edge not found in triangulated ground: "
            //           << vids.front() << ", " << vids.back() << std::endl;
            continue;
        }
        m_new.edge_data(eid_tri).label = m_old.edge_data(eid).label;

        // if the edge corresponds to a street, set the z-coordinate of the vertices
        // perchè la z delle strade è sbagliata?
        if (m_new.edge_data(eid_tri).label != -1) {
            for (uint vid : m_new.edge_vert_ids(eid_tri)) {
                m_new.vert(vid).z() = get_avg_z(m_new, vid);
            }
        }
    }
}

#endif // STREETS
