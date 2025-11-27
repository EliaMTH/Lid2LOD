#ifndef DRAW_MESH_H
#define DRAW_MESH_H

#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/gl/volume_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_sphere.h>

using namespace cinolib;

template<class M, class V, class E, class P> inline
void draw_mesh(const AbstractPolygonMesh<M,V,E,P> &m, const bool close = true)
{
    DrawablePolygonmesh<M,V,E,P> dm(m.vector_verts(), m.vector_polys());
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        dm.vert_data(vid) = m.vert_data(vid);
    }
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        dm.poly_data(pid) = m.poly_data(pid);
    }
    dm.updateGL();

    GLcanvas gui(1000, 1000);
    gui.push(&dm);

    // DrawableSphere s(m.vert(m.num_verts()-1), 1, Color::BLUE());
    // gui.push(&s);

    SurfaceMeshControls<DrawablePolygonmesh<M,V,E,P>> menu(&dm, &gui, "m");
    gui.push(&menu);
    gui.launch();

    if (close) exit(0);
}

template<class M, class V, class E, class F, class P> inline
void draw_mesh(const AbstractPolyhedralMesh<M,V,E,F,P> &m, const bool close = true)
{
    DrawablePolyhedralmesh<M,V,E,F,P> dm(m.vector_verts(), m.vector_faces(), m.vector_polys(), m.vector_windings());
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        dm.vert_data(vid) = m.vert_data(vid);
    }
    for (uint eid=0; eid<m.num_edges(); ++eid) {
        dm.edge_data(eid) = m.edge_data(eid);
    }
    for (uint fid=0; fid<m.num_faces(); ++fid) {
        dm.face_data(fid) = m.face_data(fid);
    }
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        dm.poly_data(pid) = m.poly_data(pid);
    }
    dm.updateGL();

    GLcanvas gui(1000, 1000);
    gui.push(&dm);

    VolumeMeshControls<DrawablePolyhedralmesh<M,V,E,F,P>> menu(&dm, &gui, "m");
    gui.push(&menu);
    gui.launch();

    if (close) exit(0);
}

#endif // DRAW_MESH_H
