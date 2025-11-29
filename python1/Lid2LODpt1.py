import os
import sys
import numpy as np
import shapefile
import laspy
from matplotlib.path import Path  # for inpolygon equivalent
from scipy.spatial import cKDTree

## -------------------------------------------------------------------------------
## -------------------- Subfunctions (Helper/Utility Functions) ------------------
## -------------------------------------------------------------------------------

def read_las_file(file_path):
    """
    Reads a LAS file and extracts XYZ coordinates and classification.

    Parameters:
        file_path (str): Path to the .las file.

    Returns:
        coords (np.ndarray): Nx3 array of XYZ coordinates.
        classification (np.ndarray): N-length array of classifications.
    """
    # Open the LAS file
    las = laspy.read(file_path)
    
    # Extract XYZ coordinates as Nx3 array
    coords = np.vstack((las.x, las.y, las.z)).T
    
    # Extract classification
    classification = las.classification

    classification = np.array(las.classification)
    
    return coords, classification


##

def read_shapefile_to_dict(shapefile_path):
    """
    Reads a shapefile and returns a MATLAB-style structure array with NaNs
    separating parts of polygons/polylines, similar to MATLAB's shaperead.
    
    Parameters:
        shp_path (str): Path to the shapefile (without .shp extension required)
        
    Returns:
        M (list of dict): Each element corresponds to a shape record with:
            - 'Geometry': shape type ('Polygon', 'Polyline', etc.)
            - 'X': x-coordinates with NaN separators
            - 'Y': y-coordinates with NaN separators
            - other attributes from the shapefile record
    """
    sf = shapefile.Reader(shapefile_path)
    M = []
    
    for shape_rec in sf.shapeRecords():
        shape = shape_rec.shape
        rec = shape_rec.record.as_dict() if hasattr(shape_rec.record, "as_dict") else dict(zip(sf.fields[1:], shape_rec.record))
        
        # Initialize output coordinates
        x_coords = []
        y_coords = []
        
        # Handle only shapes with points (ignore NullShapes)
        if shape.points:
            parts = list(shape.parts) + [len(shape.points)]  # add end index
            for i in range(len(shape.parts)):
                start_idx = parts[i]
                end_idx = parts[i+1]
                part = shape.points[start_idx:end_idx]
                if part:
                    xs, ys = zip(*part)
                    x_coords.extend(xs)
                    y_coords.extend(ys)
                    # Insert NaN between parts
                    x_coords.append(np.nan)
                    y_coords.append(np.nan)
        
        # Build MATLAB-style structure
        shape_struct = {
            'Geometry': shape.shapeTypeName,  # e.g., 'Polygon', 'Polyline'
            'X': np.array(x_coords),
            'Y': np.array(y_coords)
        }
        
        # Add attributes
        for k, v in rec.items():
            shape_struct[k] = v
        
        M.append(shape_struct)
    
    return M

def ReduceScopeWhatDoIDo1(xyz = None,M_i = None,other = None): 
    # building and numpy arrays
    X = np.array(M_i['X'], dtype=float)
    Y = np.array(M_i['Y'], dtype=float)


    # ----
    # use_these starts as all True
    use_these = np.ones_like(X, dtype=bool)

    # Loop from index 1 to end (since MATLAB is 1-based)
    for j in range(1, len(X)):
        if not np.isnan(X[j]) and X[j] == X[j-1] and Y[j] == Y[j-1]:
            use_these[j] = False

    X = X[use_these]
    Y = Y[use_these]
    # ----

    control = np.where(np.isnan(X))[0]
    xy = np.column_stack((X, Y))  # shape will be (n, 2)
    xMin = np.nanmin(X)
    xMax = np.nanmax(X)
    yMin = np.nanmin(Y)
    yMax = np.nanmax(Y)
    sogliaX = (xMax - xMin) / 500
    sogliaY = (yMax - yMin) / 500

    # Taglio intorno di bbx per semplificare il rangesearch

    mask = (
    (xyz[:, 0] < xMax + sogliaX) & (xyz[:, 0] > xMin - sogliaX) &
    (xyz[:, 1] < yMax + sogliaY) & (xyz[:, 1] > yMin - sogliaY)
    )

    building = xyz[mask, :]

    sogliaX_aus = (xMax - xMin) / 4
    sogliaY_aus = (yMax - yMin) / 4
    mask_aus = (
    (other[:, 0] < xMax + sogliaX_aus) & (other[:, 0] > xMin - sogliaX_aus) &
    (other[:, 1] < yMax + sogliaY_aus) & (other[:, 1] > yMin - sogliaY_aus)
    )

    other_aus = other[mask_aus, :]

    # Initialize label array
    label = np.full(X.size, 1, dtype=float)  # i*ones(numel(X),2)
    lab = 0

    # Fill second column of label
    for k in range(X.size):  # Python indexing starts at 0
        if np.any(control == k): 
            lab += 1
        label[k] = lab

    # Select xy points where label's second column is 0
    pip = xy[label == 0, :]

    if pip.shape[0] > 1:  # Need at least 2 points for inpolygon
        # Create a polygon path
        poly_path = Path(pip[:-1, :])  # MATLAB pip(1:end-1,:) → pip[:-1,:] in Python
        # Check which points in building are inside the polygon
        in_mask = poly_path.contains_points(building[:, :2])
        if np.any(in_mask):
            building = building[in_mask, :]


    return other_aus,X,Y,building

##

def DefineLabels(X, i):
    """
    Args:
        X (array-like): Input array of numbers (can contain NaN)
        i (int/float): Label identifier
    
    Returns:
        np.ndarray: lab array with shape (N, 2)
    """
    X = np.array(X, dtype=float)
    lab = []
    l = 0

    for x in X:
        if np.isnan(x):
            if lab:  # only remove last row if lab is not empty
                lab.pop()
            l += 1
        else:
            lab.append([i, l])
    
    lab = np.array(lab)
    return lab

##

def DefineTR(lab, N_abs):
    """   
    Args:
        lab (np.ndarray): lab array of shape (N, 2)
        N_abs (int): offset to be added to indices
    
    Returns:
        np.ndarray: tr array (Nx4) with 0-based indexing
    """
    tr = []
    N = 0

    for k in range(lab[:, 1].max() + 1):  # loop over unique labels
        indices_k = np.where(lab[:, 1] == k)[0]
        N_aus = len(indices_k) - 1

        for s in range(1, N_aus+1):  # MATLAB 1:N_aus-1
            tr.append([N + s + N_abs, N + s + N_abs + 1, s + 1 + N, s + N])  # careful with offsets
            # Python is 0-based, adjust if needed

        # Append the extra row after the loop
        tr.append([N + 1 + N_abs, 1 + N, N + s + 1, N+s+1+N_abs])

        N += N_aus + 1

    tr = np.array(tr, dtype=int) - 1  # convert to 0-based indexing like MATLAB
    return tr

##

def exportFacades(tr, points, dir_name):
    """
    Export mesh facades to an OFF file.

    Parameters:
    - tr: numpy array of shape (n_faces, 4) containing face indices
    - points: numpy array of shape (n_vertices, 3) containing vertex coordinates
    - dir_name: string, directory path + filename prefix
    """
    
    # Create the array like in MATLAB: first column 4, rest are tr
    aus = np.ones((tr.shape[0], 5), dtype=int) * 4
    aus[:, 1:5] = tr
    
    name_file = f"{dir_name}/facades.off"
    
    with open(name_file, 'w') as f:
        # Write OFF header
        f.write("OFF\n")
        # Write number of vertices, faces, edges
        f.write(f"{points.shape[0]} {tr.shape[0]} 0\n")
        # Write vertex coordinates
        for p in points:
            f.write(f"{p[0]:12.8f} {p[1]:12.8f} {p[2]:12.8f}\n")
        # Write faces
        for face in aus:
            f.write(f"{face[0]:5d} {face[1]:5d} {face[2]:5d} {face[3]:5d} {face[4]:5d}\n")

##

def exportPavement(points_build_f, lab, dirNames):
    aus = np.arange(points_build_f.shape[0])
    nameFile = f"{dirNames}/pavement_polygon.off"

    with open(nameFile, 'w') as fileID:
        fileID.write('OFF\n')
        fileID.write(f"{points_build_f.shape[0]} {int(np.max(lab[:, 1])) + 1} 0\n")

        # Write vertex coordinates
        for point in points_build_f:
            fileID.write(f"{point[0]:12.8f} {point[1]:12.8f} {point[2]:12.8f}\n")

        # Write face indices
        max_label = int(np.max(lab[:, 1]))
        for k in range(max_label + 1):
            indices = aus[lab[:, 1] == k]
            fileID.write(f"{len(indices)} ")
            
            if k > 0:
                for idx in indices:
                    fileID.write(f"{idx} ")
            else:
                for idx in indices[::-1]:  # reverse for k=0
                    fileID.write(f"{idx} ")
            
            fileID.write('\n')
        fileID.write('\n')

##

def exportRoof(points_build_r, lab, dirNames):
    aus = np.arange(points_build_r.shape[0])
    nameFile = f"{dirNames}/roof_polygon.off"
    
    with open(nameFile, 'w') as fileID:
        # Write the OFF header
        fileID.write('OFF\n')
        
        # Write counts: number of vertices, number of faces, 0
        fileID.write(f"{points_build_r.shape[0]} {int(np.max(lab[:,1]))+1} 0\n")
        
        # Write the vertex coordinates
        for point in points_build_r:
            fileID.write(f"{point[0]:.8f} {point[1]:.8f} {point[2]:.8f}\n")
        
        # Write the face indices
        for k in range(int(np.max(lab[:,1]))+1):
            indices = aus[lab[:,1] == k]
            fileID.write(f"{len(indices)} ")
            
            if k > 0:
                for idx in indices:
                    fileID.write(f"{idx} ")
            else:
                for idx in indices[::-1]:
                    fileID.write(f"{idx} ")
            fileID.write('\n')

## -------------------------------------------------------------------------------
## ---------------------------------- MAIN -------------------------------
## -------------------------------------------------------------------------------

def main(building_footprints,las_path,temp_fold):

    ## Read pointcloud and assign attributes

    xyz, pt_classification = read_las_file(las_path) # Should already be numpy array

    other = xyz[pt_classification == 2, :].astype(float)
    xyz = xyz[pt_classification == 6, :].astype(float)


    ## Read shapefile
    M = read_shapefile_to_dict(building_footprints)
    ## Initialize further variables
    zMin_default = np.mean(other[:, 2])

    for i in range(0, len(M)):

        # Computing other_aus, control, X, Y, xy, zMin, building
        other_aus,X,Y,building = ReduceScopeWhatDoIDo1(xyz,M[i],other)

        if building.shape[0] > X.shape[0]:
            if building.shape[0] > 50000:
                building = building[::10, :]

            xy = []
            for x, y in zip(X, Y):
                if np.isnan(x):
                    if xy:  # only remove last row if xy is not empty
                        xy.pop()
                else:
                    xy.append([x, y])
            
            xy = np.array(xy)  # convert list to numpy array


            # Step 1 Organize Data
            lab = DefineLabels(X,i)

            if other_aus.shape[0] > 0:
                tree = cKDTree(other_aus[:, :2])
                _, I = tree.query(xy)
                zMin = np.min(other_aus[I, 2])
            else:
                zMin = zMin_default
            zMax = np.max(building[:, 2])




            points_aus_ground = np.zeros((xy.shape[0],6))

            z = np.full((xy.shape[0], 1), zMin)
            xyz_aus = np.hstack([xy, z])
            points_aus_ground[:, 0:3] = xyz_aus
            points_aus_ground[:, 4] = i
            
            z = np.full((xy.shape[0], 1), zMax)
            points_aus_roof = np.ones((xy.shape[0], 6))
            xyz_aus = np.hstack([xy, z])
            points_aus_roof[:, 0:3] = xyz_aus
            points_aus_roof[:, 4] = i
            
            points_build_f = points_aus_ground[:, 0:3]
            points_build_r = points_aus_roof[:, 0:3]
            points = np.vstack([points_build_f, points_build_r])
            
            tr = DefineTR(lab, points_build_f.shape[0])

            # Export OFF files
            dirNames = f'{temp_fold}/building{i+1}'
            os.makedirs(dirNames, exist_ok=True)
            exportFacades(tr,points,dirNames)
            exportPavement(points_build_f,lab,dirNames)
            exportRoof(points_build_r,lab,dirNames)
    
    print("Saved building polygons:", len(M))

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python Lid2LODpt1.py <building_footprints.shp> <points.las> <temp_output_folder>")
        exit(1)

    building_footprints = sys.argv[1]
    las_path = sys.argv[2]
    temp_fold = sys.argv[3]

    main(building_footprints, las_path, temp_fold)
