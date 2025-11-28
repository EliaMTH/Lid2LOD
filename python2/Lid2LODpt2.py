import os
import sys
import json
import numpy as np


## -------------------------------------------------------------------------------
## -------------------- Subfunctions (Helper/Utility Functions) ------------------
## -------------------------------------------------------------------------------

def off_to_dict(filepath):
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    assert lines[0] == "OFF", "Not a valid OFF file"
    
    # counts
    n_vertices, n_faces, _ = map(int, lines[1].split())
    
    # vertices
    vertices = []
    for i in range(2, 2 + n_vertices):
        x, y, z = map(float, lines[i].split())
        vertices.append([x, y, z])
    
    # faces (remove the first entry, which is just count)
    faces = []
    for i in range(2 + n_vertices, 2 + n_vertices + n_faces):
        parts = list(map(int, lines[i].split()))
        indices = parts[1:]  # skip the first value (number of vertices)
        faces.append(indices)
    
    # geometry formatting
    if len(faces) == 1:
        geometry = [faces[0]]   # single polygon
    else:
        geometry = [faces]      # multiple polygons
    
    result = {
        "geometry": geometry,
        "vertices": vertices
    }
    
    return result


def merge_off_models(models, merged_name="merged_buildings"):
    all_vertices = []
    all_geometry = []

    vertex_offset = 0
    for model in models:
        for g in model["geometry"]:
            # g could be a polygon (list of ints) or a multi-ring (list of list of ints)
            if all(isinstance(el, int) for el in g):
                # single polygon
                geom_shifted = [i + vertex_offset for i in g]
                all_geometry.append(geom_shifted)
            else:
                # multi-ring polygon, shift each ring separately
                geom_shifted = [[i + vertex_offset for i in ring] for ring in g]
                all_geometry.append(geom_shifted)
        all_vertices.extend(model["vertices"])
        vertex_offset += len(model["vertices"])

    # Deduplicate vertices
    unique_vertices = []
    vertex_map = {}
    for i, v in enumerate(all_vertices):
        v_tuple = tuple(round(c, 8) for c in v)
        if v_tuple not in vertex_map:
            vertex_map[v_tuple] = len(unique_vertices)
            unique_vertices.append(list(v_tuple))

    # Update geometry indices
    new_geometry = []
    for g in all_geometry:
        if all(isinstance(el, int) for el in g):
            new_geometry.append([vertex_map[tuple(round(all_vertices[i][j], 8) for j in range(3))] for i in g])
        else:
            new_geometry.append([
                [vertex_map[tuple(round(all_vertices[i][j], 8) for j in range(3))] for i in ring]
                for ring in g
            ])

    merged_model = {
        "id": merged_name,
        "lod": "1.1",
        "geometry": new_geometry,
        "vertices": unique_vertices
    }

    return merged_model

import os

def off_to_dict_walls(filepath):
    with open(filepath, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    assert lines[0] == "OFF", "Not a valid OFF file"
    
    # counts
    n_vertices, n_faces, _ = map(int, lines[1].split())
    
    # vertices
    vertices = []
    for i in range(2, 2 + n_vertices):
        x, y, z = map(float, lines[i].split())
        vertices.append([x, y, z])
    
    # faces (remove the first entry, which is just count)
    faces = []
    for i in range(2 + n_vertices, 2 + n_vertices + n_faces):
        parts = list(map(int, lines[i].split()))
        indices = parts[1:]  # skip the first number
        faces.append(indices)
    
    # geometry = keep all faces separate, no nesting
    geometry = faces
    
    result = {
        "geometry": geometry,
        "vertices": vertices
    }
    
    return result

def parse_OFF_building(building_path):
    buildings_data_f = off_to_dict(f"{building_path}/pavement_polygon.off")
    buildings_data_r = off_to_dict(f"{building_path}/roof_polygon.off")
    buildings_data_w = off_to_dict_walls(f"{building_path}/facades.off")
    buildings_data = merge_off_models([buildings_data_f, buildings_data_r, buildings_data_w], os.path.basename(building_path))

    return buildings_data


def parse_OFF_building_dataset(main_folder):
    """
    Parses all subfolders in a main folder recursively, 
    running parse_OFF_building on folders containing the required OFF files.

    Args:
        main_folder (str): Path to the main folder containing building subfolders.

    Returns:
        dict: Dictionary mapping subfolder names to their parsed building data.
    """
    all_buildings_data = {}

    for root, dirs, files in os.walk(main_folder):
        # Check for required files
        pavement_file = "pavement_polygon.off"
        facades_file = "facades.off"
        building_pavement_file = "pavement_polygon.off"

        if (pavement_file in files and
            os.path.exists(os.path.join(root, building_pavement_file)) and
            os.path.exists(os.path.join(root, facades_file))):
            try:
                data = parse_OFF_building(root)
                all_buildings_data[os.path.basename(root)] = data
            except Exception as e:
                print(f"Failed to parse {root}: {e}")

    return all_buildings_data

def create_cityjson(buildings=None, output_file="city.json"):

    # Global vertex list and mapping to deduplicate
    global_vertices = []
    vertex_map = {}

    # def get_vertex_index(vertex):
    #     vtuple = tuple(vertex)
    #     if vtuple not in vertex_map:
    #         vertex_map[vtuple] = len(global_vertices)
    #         global_vertices.append(vertex)
    #     return vertex_map[vtuple]
    
    def get_vertex_index(vertex, vertices):
        """Map a vertex or index to a global index."""
        if isinstance(vertex, int):
            vertex = vertices[vertex]  # get the actual coordinate
        vtuple = tuple(vertex)
        if vtuple not in vertex_map:
            vertex_map[vtuple] = len(global_vertices)
            global_vertices.append(vertex)
        return vertex_map[vtuple]

    # Base CityJSON structure
    cityjson_dict = {
        "type": "CityJSON",
        "version": "1.1",
        "CityObjects": {},
        "vertices": [],
    }

    # Build CityObjects
    for b in buildings:
        geom_faces = []
        semantic_surfaces = []
        semantics_values = []

        for idx, face in enumerate(b["geometry"]):
            if all(isinstance(f, list) for f in face):
                # Multiple rings (e.g., floor with holes)
                face_indices = [[get_vertex_index(v, b["vertices"]) for v in ring] for ring in face]
            else:
                # Single ring
                face_indices = [[get_vertex_index(v, b["vertices"]) for v in face]]


            geom_faces.append({"type": "Polygon", "boundaries": face_indices})

            # Determine semantic type
            if idx == 0:
                semantic_type = "FloorSurface"
            elif idx == 1:
                semantic_type = "RoofSurface"
            else:
                semantic_type = "WallSurface"

            semantic_surfaces.append({"type": semantic_type})
            semantics_values.append(idx)

        cityjson_dict["CityObjects"][b["id"]] = {
            "type": "Building",
            "geometry": [{
                "type": "Solid",
                "lod": b.get("lod", "1.2"),
                "boundaries": [[face["boundaries"] for face in geom_faces]],
                "semantics": {
                    "surfaces": semantic_surfaces,
                    "values": [semantics_values]
                }
            }]
        }

    # Add vertices
    cityjson_dict["vertices"] = global_vertices

    # Compute transform (required by CityJSON 1.1)
    vertices_array = np.array(global_vertices)
    min_coords = vertices_array.min(axis=0)
    scale = [0.01, 0.01, 0.01]
    translate = min_coords.tolist()

    cityjson_dict["transform"] = {
        "scale": scale,
        "translate": translate
    }

    # Save as file
    if not output_file.endswith(".json"):
        output_file = f"{output_file}.city.json"
    else:
        output_file = output_file[:-5] + ".city.json"

    with open(output_file, "w") as f:
        json.dump(cityjson_dict, f, indent=2)

    print(f"CityJSON written to {output_file}")

## -------------------------------------------------------------------------------
## ---------------------------------- MAIN -------------------------------
## -------------------------------------------------------------------------------

def main(working_folder,outpath):
    print
    bb = parse_OFF_building_dataset(working_folder)
    create_cityjson(list(bb.values()), outpath)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 Lid2LODpt2.py <working_folder> <outpath>")
        sys.exit(1)

    working_folder = sys.argv[1]
    outpath = sys.argv[2]

    main(working_folder, outpath)