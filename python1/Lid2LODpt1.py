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


    return other_aus,control,X,Y,xy,building

##

def SelectBuildingsAndLabels(building, i, X, xy, control):
    # Ensure building is float type
    building = building.astype(float)

    # Initialize label array
    label = np.full((X.size, 2), i, dtype=float)  # i*ones(numel(X),2)
    lab = 0

    # Fill second column of label
    for k in range(X.size):  # Python indexing starts at 0
        if np.any(control == k): 
            lab += 1
        label[k, 1] = lab

    # Select xy points where label's second column is 0
    pip = xy[label[:, 1] == 0, :]

    if pip.shape[0] > 1:  # Need at least 2 points for inpolygon
        # Create a polygon path
        poly_path = Path(pip[:-1, :])  # MATLAB pip(1:end-1,:) → pip[:-1,:] in Python
        # Check which points in building are inside the polygon
        in_mask = poly_path.contains_points(building[:, :2])
        if np.any(in_mask):
            building = building[in_mask, :]

    return building

##

def Define2DProjections(M_i):
    X = np.array(M_i['X'], dtype=float)
    Y = np.array(M_i['Y'], dtype=float)
    
    xy = []
    
    for x, y in zip(X, Y):
        if np.isnan(x):
            if xy:  # only remove last row if xy is not empty
                xy.pop()
        else:
            xy.append([x, y])
    
    xy = np.array(xy)  # convert list to numpy array
    
    return X, xy

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

        N += N_aus

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

## Elaborate Polygons
def ElaboratePoligons(k, xy, control, num_piano, building, other_aus, zMin, 
                       points_LOD2_facade_ground, points_LOD2_facade_roof):
    # Check if k+1 is in control
    if not np.any(k + 1 == control):
        n_rows = xy.shape[0]

        # Determine XY1 and XY2
        if k == n_rows - 1 and np.isscalar(control):
            XY1 = xy[k, :]
            XY2 = xy[0, :]
        elif k == n_rows - 1 and len(control) > 1:
            XY1 = xy[k, :]
            XY2 = xy[control[-2] + 1, :]
        else:
            XY1 = xy[k, :]
            XY2 = xy[k + 1, :]

        # Create linearly spaced coordinates
        if abs(XY1[0] - XY2[0]) < 0.1:
            y_piano = np.linspace(min(XY1[1], XY2[1]), max(XY1[1], XY2[1]), num_piano)
            x_piano = np.full(num_piano, XY1[0])
        else:
            if abs(XY1[1] - XY2[1]) < 0.1:
                x_piano = np.linspace(min(XY1[0], XY2[0]), max(XY1[0], XY2[0]), num_piano)
                y_piano = np.full(num_piano, XY1[1])
            else:
                m = (XY2[1] - XY1[1]) / (XY2[0] - XY1[0])
                if not np.isnan(m):
                    q = XY1[1] - m * XY1[0]
                    x_piano = np.linspace(min(XY1[0], XY2[0]), max(XY1[0], XY2[0]), num_piano)
                    y_piano = m * x_piano + q
                else:
                    return points_LOD2_facade_ground, points_LOD2_facade_roof  # exit function if plane scenario is not accepted

        # Initialize piano points
        piano = np.zeros((len(x_piano) * len(y_piano), 3))

        # Build a KD-tree for nearest neighbor search
        if other_aus.shape[0] > 0:
            tree = cKDTree(other_aus[:, :2])

        for s in range(len(x_piano)):
            idx_start = s * len(y_piano)
            idx_end = (s + 1) * len(y_piano)
            piano[idx_start:idx_end, 0] = x_piano[s]
            piano[idx_start:idx_end, 1] = y_piano

            # Find nearest neighbor in other_aus
            if other_aus.shape[0] > 0:
                _, I = tree.query([x_piano[s], y_piano[s]])
                zMin = np.min(other_aus[I, 2])

            zMax = np.max(building[:, 2])
            zAus = np.linspace(zMin, zMax, num_piano)
            piano[idx_start:idx_end, 2] = zAus[::-1]

            # Update facade points
            if s == 0 or s == len(x_piano) - 1:
                if points_LOD2_facade_ground.shape[0] > 0:
                    # Check if point already exists
                    dists = np.linalg.norm(points_LOD2_facade_roof[:, :2] - np.array([x_piano[s], y_piano[s]]), axis=1)
                    if np.any(dists < 1e-5):
                        index = np.where((points_LOD2_facade_ground[:, 0] == x_piano[s]) &
                                         (points_LOD2_facade_ground[:, 1] == y_piano[s]))[0]
                        if len(index) > 0 and points_LOD2_facade_roof[index, 2] > zMax:
                            points_LOD2_facade_roof[index, 2] = zMax
                    else:
                        points_LOD2_facade_ground = np.vstack([points_LOD2_facade_ground, [x_piano[s], y_piano[s], zMin]])
                        points_LOD2_facade_roof = np.vstack([points_LOD2_facade_roof, [x_piano[s], y_piano[s], zMax]])
                else:
                    points_LOD2_facade_ground = np.vstack([points_LOD2_facade_ground, [x_piano[s], y_piano[s], zMin]])
                    points_LOD2_facade_roof = np.vstack([points_LOD2_facade_roof, [x_piano[s], y_piano[s], zMax]])

    return points_LOD2_facade_ground, points_LOD2_facade_roof

## -------------------------------------------------------------------------------
## ---------------------------------- MAIN -------------------------------
## -------------------------------------------------------------------------------

def main(building_footprints,las_path,temp_fold):

    ## Read pointcloud and assign attributes

    xyz, pt_classification = read_las_file(las_path) # Should already be numpy array

    other = xyz[pt_classification == 2, :].astype(float)
    # veget = xyz[pt_classification == 4, :].astype(float)
    xyz = xyz[pt_classification == 6, :].astype(float)


    ## Read shapefile
    M = read_shapefile_to_dict(building_footprints)
    ## Initialize further variables
    num_piano = 10
    zMin = np.mean(other[:, 2])

    for i in range(0, len(M)):
        # Computing other_aus, control, X, Y, xy, zMin, building
        other_aus,control,X,Y,xy,building = ReduceScopeWhatDoIDo1(xyz,M[i],other)

        if building.shape[0] > 60:
            if building.shape[0] > 50000:
                building = building[::10, :]
            if building.shape[0] > 0:
                # Compute points_LOD2_*
                building = SelectBuildingsAndLabels(building,i,X,xy,control)
                points_LOD2_facade_ground = np.empty((0, 3))
                points_LOD2_facade_roof = np.empty((0, 3))
                for k in range(xy.shape[0] - 1):
                    points_LOD2_facade_ground,points_LOD2_facade_roof = ElaboratePoligons(k,xy,control,num_piano,building,other_aus,zMin,points_LOD2_facade_ground,points_LOD2_facade_roof)

                # Step 1 Organize Data
                X,xy = Define2DProjections(M[i])
                lab = DefineLabels(X,i)
                points_aus_ground = np.zeros((xy.shape[0],6))

                z = np.full((xy.shape[0], 1), np.mean(points_LOD2_facade_ground[:, 2]))
                xyz_aus = np.hstack([xy, z])
                points_aus_ground[:, 0:3] = xyz_aus
                points_aus_ground[:, 4] = i
                
                z = np.full((xy.shape[0], 1), np.max(points_LOD2_facade_roof[:, 2]))
                points_aus_roof = np.ones((xy.shape[0], 6))
                xyz_aus = np.hstack([xy, z])
                points_aus_roof[:, 0:3] = xyz_aus
                points_aus_roof[:, 4] = i
                
                points_build_f = points_aus_ground[:, 0:3]
                points_build_r = points_aus_roof[:, 0:3]
                points = np.vstack([points_build_f, points_build_r])
                
                tr = DefineTR(lab, points_build_f.shape[0])

                # Export OFF files
                dirNames = f'./{temp_fold}/building{i+1}'
                os.makedirs(dirNames, exist_ok=True)
                exportFacades(tr,points,dirNames)
                exportPavement(points_build_f,lab,dirNames)
                exportRoof(points_build_r,lab,dirNames)

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python Lid2LODpt1.py <building_footprints.shp> <points.las> <temp_output_folder>")
        exit(1)

    building_footprints = sys.argv[1]
    las_path = sys.argv[2]
    temp_fold = sys.argv[3]

    main(building_footprints, las_path, temp_fold)
