#ifndef MERGE_MESH_AT_COINCIDENT_VERTICES_TOLL_H
#define MERGE_MESH_AT_COINCIDENT_VERTICES_TOLL_H

#include <cinolib/meshes/meshes.h>
#include <cinolib/octree.h>

using namespace cinolib;

template<class M, class V, class E, class P>
CINO_INLINE void merge_meshes_at_coincident_vertices(const AbstractPolygonMesh<M,V,E,P> & m1,
                                                     const AbstractPolygonMesh<M,V,E,P> & m2,
                                                     AbstractPolygonMesh<M,V,E,P> & res,
                                                     const double proximity_thresh)
{
    Octree octree;
    for(uint vid=0; vid<m1.num_verts(); ++vid)
        octree.push_sphere(vid, m1.vert(vid), proximity_thresh);

    res = m1;

    std::map<uint,uint> vmap;
    for(uint vid=0; vid<m2.num_verts(); ++vid)
    {
        vec3d p = m2.vert(vid);

        std::unordered_set<uint> ids;
        if(octree.contains(p, false, ids))
        {
            // WARNING: I am assuming that the mapping is one to one at most
            assert(ids.size()==1);
            vmap[vid] = *ids.begin();
        }
        else
        {
            uint fresh_id = res.vert_add(p);
            vmap[vid] = fresh_id;
        }
    }

    for(uint pid=0; pid<m2.num_polys(); ++pid)
    {
        auto p = m2.poly_verts_id(pid);

        for(auto & vid : p)
            vid = vmap.at(vid);

        int test_id = res.poly_id(p);
        if(test_id<0)
        {
            uint fresh_id = res.poly_add(p);
        }
    }
}

#endif // MERGE_MESH_AT_COINCIDENT_VERTICES_TOLL_H
