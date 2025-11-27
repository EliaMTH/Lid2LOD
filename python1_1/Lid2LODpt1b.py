import numpy as np
import os
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull
import laspy

def minBoundingBox(X):
    """
    Fully faithful translation of the MATLAB minBoundingBox.m
    Input:  X is a 2×n numpy array
    Output: bb is a 2×4 numpy array of bounding box corners
    """

    # ---- convex hull ----
    points = X.T  # SciPy ConvexHull expects Nx2
    hull = ConvexHull(points)
    CH = X[:, hull.vertices]   # 2 × k

    # ---- angles of convex hull edges ----
    E = np.diff(CH, axis=1)    # 2 × (k-1)
    T = np.arctan2(E[1, :], E[0, :])
    T = np.mod(T, np.pi/2)
    T = np.unique(T)

    # ---- construct the R matrix exactly as MATLAB does ----
    # MATLAB code:
    # R = cos( reshape(repmat(T,2,2),2*length(T),2)
    #          + repmat([0 -pi ; pi 0]/2,length(T),1) );
    #
    # Build A = repmat(T, 2, 2) then reshape
    A = np.tile(T, (2, 2))  # shape (2, 2*len(T))
    A = np.reshape(A, (2*len(T), 2), order='F')

    # B = repmat([0 -pi; pi 0]/2, length(T), 1)
    B = np.array([[0, -np.pi], [np.pi, 0]]) / 2
    B = np.tile(B, (len(T), 1))

    R = np.cos(A + B)  # shape (2*len(T), 2)

    # ---- rotate convex hull CH using all angles ----
    RCH = R @ CH  # (2*len(T)) × k
    RCH = RCH

    # ---- compute bounding sizes for each rotation ----
    bsize = np.max(RCH, axis=1) - np.min(RCH, axis=1)  # length 2*nAngles
    bsize = bsize.reshape(2, -1, order='F')            # 2 × nAngles
    area = np.prod(bsize, axis=0)                      # nAngles

    # ---- find minimal area ----
    i = np.argmin(area)

    # ---- compute bounding box in rotated frame ----
    # Rf = R(2*i+[-1 0], :)
    row1 = 2*i - 2 
    row2 = 2*i - 1
    Rf = np.vstack((R[row1, :], R[row2, :]))           # 2×2

    bound = Rf @ CH
    bmin = np.min(bound, axis=1)
    bmax = np.max(bound, axis=1)

    Rf = Rf.T

    # ---- compute final bounding box corners ----
    # MATLAB:
    # bb(:,4) = bmax(1)*Rf(:,1) + bmin(2)*Rf(:,2)
    # bb(:,1) = bmin(1)*Rf(:,1) + bmin(2)*Rf(:,2)
    # bb(:,2) = bmin(1)*Rf(:,1) + bmax(2)*Rf(:,2)
    # bb(:,3) = bmax(1)*Rf(:,1) + bmax(2)*Rf(:,2)
    #
    bb = np.zeros((2, 4))

    bb[:, 3] = bmax[0] * Rf[:, 0] + bmin[1] * Rf[:, 1]
    bb[:, 0] = bmin[0] * Rf[:, 0] + bmin[1] * Rf[:, 1]
    bb[:, 1] = bmin[0] * Rf[:, 0] + bmax[1] * Rf[:, 1]
    bb[:, 2] = bmax[0] * Rf[:, 0] + bmax[1] * Rf[:, 1]

    return bb

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


# --------------------------------------------------
# ------------------- MAIN ---------------------
# --------------------------------------------------

def main(las_path,outname):
    xyz, pt_classification = read_las_file(las_path)

    other = xyz[pt_classification == 2].astype(float)
    xyz   = xyz[pt_classification == 6].astype(float)

    #------------------------------------------
    # Compute min bounding box on (other; xyz)
    #------------------------------------------
    aus = np.vstack((other, xyz))
    mBB = minBoundingBox(aus[:, :2].T)      # 4×2 array like MATLAB mBB'

    # Rearrange to match MATLAB structure  
    # xy = np.zeros((mBB.shape[1], 3))
    # xy[:, :2] = mBB
    xy = np.hstack((mBB.T, np.zeros((mBB.shape[1], 1))))


    #------------------------------------------
    # Expand corners by T = 100
    #------------------------------------------
    T = 100
    xy_mod = xy.copy()

    xy_mod[0,0] += T;  xy_mod[0,1] -= T
    xy_mod[1,0] -= T;  xy_mod[1,1] -= T
    xy_mod[2,0] -= T;  xy_mod[2,1] += T
    xy_mod[3,0] += T;  xy_mod[3,1] += T

    xy = xy_mod

    #------------------------------------------
    # Sample the boundary (num_piano = 15)
    #------------------------------------------
    num_piano = 15
    bound = []

    for s in range(xy.shape[0] - 1):
        x1, y1 = xy[s, :2]
        x2, y2 = xy[s+1, :2]

        if abs(x1 - x2) > 1.5:
            x = np.linspace(min(x1, x2), max(x1, x2), num_piano)
            y = np.full(num_piano, y1)
        else:
            y = np.linspace(min(y1, y2), max(y1, y2), num_piano)
            x = np.full(num_piano, x1)

        bound.append(np.column_stack((x, y)))

    # Last segment
    x1, y1 = xy[-1, :2]
    x2, y2 = xy[0, :2]

    if abs(x1 - x2) > 1.5:
        x = np.linspace(min(x1, x2), max(x1, x2), num_piano)
        y = np.full(num_piano, y1)
    else:
        y = np.linspace(min(y1, y2), max(y1, y2), num_piano)
        x = np.full(num_piano, x1)

    bound.append(np.column_stack((x, y)))

    bound = np.vstack(bound)
    bound = np.column_stack((bound, np.zeros(len(bound))))

    xy = bound

    #------------------------------------------
    # KNN height assignment (MATLAB knnsearch)
    #------------------------------------------
    tree = KDTree(other[:, :2])

    for i in range(xy.shape[0]):
        dist, index = tree.query(xy[i, :2])
        xy[i, 2] = other[index, 2]

    # Reverse order like MATLAB (end:-1:1)
    xy = xy[::-1, :]

    #------------------------------------------
    # Write OFF file
    #------------------------------------------
    nameFile = outname

    with open(nameFile, 'w') as f:
        f.write("OFF\n")
        f.write(f"{xy.shape[0]} 1 0\n")

        for row in xy:
            f.write(f"{row[0]:.8f} {row[1]:.8f} {row[2]:.8f}\n")

        f.write(f"{xy.shape[0]} ")
        for k in range(xy.shape[0]):
            f.write(f"{k} ")
        f.write("\n")



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute boundary OFF file from LAS input.")
    parser.add_argument("las_path", type=str, help="Path to input .las file")
    parser.add_argument("outname", type=str, help="Output OFF filename")
    args = parser.parse_args()

    # Validate input file
    if not os.path.isfile(args.las_path):
        raise FileNotFoundError(f"Input LAS file not found: {args.las_path}")

    # Ensure output directory exists
    outdir = os.path.dirname(args.outname)
    if outdir != "" and not os.path.exists(outdir):
        os.makedirs(outdir)

    main(args.las_path, args.outname)
