# [HEAVY WORK IN PROGRESS, FOR A CLEANER READ AND EXPERIENCE COME BACK SOON! Apologies in advance]

# Lid2LOD – ADDTITLE

This repository contains a fully reproducible Docker environment for running **Lid2LOD**, a 3-step processing pipeline that generates 3D building models from:

* **Building footprints** (`.shp` and eventual shapefile components .shx, .dbf, etc.)
* **Point cloud** (`.las`)

The container bundles all required Python scripts and the C++ `triangulate_city` module, so users can run the entire workflow without installing dependencies.

---

##  Features

* Runs on **Linux, Windows, and macOS** via Docker
* Includes all Python and C++ dependencies
* Simple one-command execution
* Automatically handles working directories
* Outputs LOD models into a user-defined folder


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
  (must include `.shp`, `.dbf`, `.shx`)
* A `.las` point cloud file [ADD REQUIREMENTS]
* Any optional files your scripts need

Create an output folder (or let Docker create it).


---

### **Run command (Windows)**

```powershell
docker run --rm -v "C:\path\to\my_data:/data" lid2lod_v1 `
/data/buildings.shp /data/points.las /data/output
```

---

### **Run command (Linux/macOS)**

```bash
docker run --rm -v "/path/to/my_data:/data" lid2lod_v1 \
/data/<buildings.shp> /data/<points.las> /data/<output_folder>
```

## Output

The output folder should contain:

* triangulated city output
* LOD geometries
