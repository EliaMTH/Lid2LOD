# Lid2LOD – Generating LOD1 Urban Models from Airborne LiDAR

This repository contains a fully reproducible Docker environment for running **Lid2LOD**, a processing pipeline that generates 3D building models from:

* a **Building footprints** file (`.shp` and eventual shapefile components .shx, .dbf, etc.),
* a **Point cloud** file (`.las`).

---

## Building
First, you should have docker on your machine. We tested this project using https://www.docker.com/products/docker-desktop/

From the root of the project run the following command:

```bash
docker build -t lid2lod_v1 .
```

**Warning**: if you edit `entrypoint.sh` or any shell script on Windows, it may be saved with **CRLF** line endings.
Linux containers expect **LF** only. Using CRLF causes Docker to fail with errors such as:

```
entrypoint.sh: not found
```

To fix this, configure your editor to save shell scripts with **LF** by default.


## Running 

Prepare a folder on your machine containing:
* A building footprint shapefile
  (should include `.shp`, `.dbf`, `.shx`)
* A `.las` point cloud file (point classification is required, automatic point classification is planned for next Lid2Lod version!)
* A folder where to store the outputs

You can then run Lid2LOD_v1 by running the following command, depending on your OS:

### **Run command (Windows)**
```powershell
docker run --rm -v "X:\path\to\my_data:/data" lid2lod_v1 /data/<buildings_footprints.shp> /data/<point_cloud.las> /data/<output_folder>
```

### **Run command (Linux/macOS)**

```bash
docker run --rm -v "/path/to/my_data:/data" lid2lod_v1 /data/<buildings_footprints.shp> /data/<point_cloud.las> /data/<output_folder>
```

## Output

After the process is completed, the output folder contains 4 files:
* A mesh containing all buildings (`buildings_mesh.off`)
* A mesh of the ground (`ground_mesh.off`)
* A mesh of the city (`city_mesh.off`)
* A CityJSON model (`city_JSON.city.json`)

## Reference - How to cite our work
This work was presented during STAG2025, with a paper by the title _LiD2LOD: Generating LOD1 Urban Models from Airborne LiDAR_. 
You can access the publication through the Eurographics Digital Library:

* [**DIGLIB proceedings page**](https://diglib.eg.org/items/85255f56-244f-4ee9-a38b-27f1fdf20d48)

* [**Direct PDF link**](https://diglib.eg.org/server/api/core/bitstreams/936212f9-a9aa-4b6c-80e2-f7bd8e092ac1/content)

If you use this work in your research, please cite it as follows:

```bibtex
@inproceedings{10.2312:stag.20251325,
  booktitle = {Smart Tools and Applications in Graphics - Eurographics Italian Chapter Conference},
  editor = {Comino Trinidad, Marc and Mancinelli, Claudio and Maggioli, Filippo and Romanengo, Chiara and Cabiddu, Daniela and Giorgi, Daniela},
  title = {{LiD2LOD: Generating LOD1 Urban Models from Airborne LiDAR}},
  author = {Sorgente, Tommaso and Moscoso Thompson, Elia and Romanengo, Chiara},
  year = {2025},
  publisher = {The Eurographics Association},
  ISSN = {2617-4855},
  ISBN = {978-3-03868-296-7},
  DOI = {10.2312/stag.20251325}
}
```


## Contacts
We appreciate any kind of feedback! Please reach us at:
tommaso.sorgente@cnr.it, chiara.romanengo@cnr.it, elia.moscosothompson@cnr.it
