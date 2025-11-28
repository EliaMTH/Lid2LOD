# [!!! WARNING: we are still completing this repository, for a smoother experience come back next week! Our apologies in advance. If you are interested, send us an email and we will notify you once we're done! Check out contacts below!]

# Lid2LOD – Generating LOD1 Urban Models from Airborne LiDAR

This repository contains a fully reproducible Docker environment for running **Lid2LOD**, a 3-step processing pipeline that generates 3D building models from:

* **Building footprints** (`.shp` and eventual shapefile components .shx, .dbf, etc.)
* **Point cloud** (`.las`)

---

## Building the Docker Image
First, you should have docker on your machine. We tested this project using https://www.docker.com/products/docker-desktop/

From the root of the project run the following command:

```bash
docker build -t lid2lod_v1 .
```

This creates a Docker image named **`lid2lod_v1`**.

---

## Running the Pipeline

Prepare a folder on your machine containing:

* A building footprint shapefile
  (should include `.shp`, `.dbf`, `.shx`)
* A `.las` point cloud file [ADD REQUIREMENTS]
* Any optional files your scripts need

You can then run Lid2LOD_v1 by running the following command, depending on your OS.

### **Run command (Windows)**
```powershell
docker run --rm -v "X:\path\to\my_data:/data" lid2lod_v1 /data/<buildings.shp> /data/<points.las> /data/<output_name>
```

### **Run command (Linux/macOS)**

```bash
docker run --rm -v "/path/to/my_data:/data" lid2lod_v1 /data/<buildings.shp> /data/<points.las> /data/<output_name>
```

## Output

The output folder should contain:

* triangulated city (`output_name.off`)
* CityJSON model (`output_name.city.json`)

## REFERENCE - How to city our work
https://diglib.eg.org/items/85255f56-244f-4ee9-a38b-27f1fdf20d48

## Contacts
We appreciate any kind of feedback! Please reach us at:
tommaso.sorgente@cnr.it, chiara.romanengo@cnr.it, elia.moscosothompson@cnr.it
