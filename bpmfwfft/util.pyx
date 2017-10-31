
cimport cython
import numpy as np
cimport numpy as np

cdef extern from "math.h":
    double sqrt(double)

@cython.boundscheck(False)
def cdistance(np.ndarray[np.float64_t, ndim=1] x, np.ndarray[np.float64_t, ndim=1] y):
    cdef int i, lmax
    cdef double d, tmp
    lmax = x.shape[0]
    d = 0.
    for i in range(lmax):
        tmp = x[i] - y[i]
        d += tmp*tmp
    return sqrt(d)

@cython.boundscheck(False)
def c_get_corner_crd(np.ndarray[np.int64_t, ndim=1] corner,
                     np.ndarray[np.float64_t, ndim=1] grid_x, 
                     np.ndarray[np.float64_t, ndim=1] grid_y,
                     np.ndarray[np.float64_t, ndim=1] grid_z):
    cdef:
        int i, j, k
        np.ndarray[np.float64_t, ndim=1] crd

    i, j, k = corner
    crd = np.array([grid_x[i], grid_y[j], grid_z[k]] , dtype=float)
    return crd

@cython.boundscheck(False)
def c_is_in_grid(np.ndarray[np.float64_t, ndim=1] atom_coordinate,
                 np.ndarray[np.float64_t, ndim=1] origin_crd,
                 np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd):
    
    cdef int i, lmax
    lmax = atom_coordinate.shape[0]
    for i in range(lmax):
        if (atom_coordinate[i] < origin_crd[i]) or (atom_coordinate[i] >= uper_most_corner_crd[i]):
            return False
    return True

@cython.boundscheck(False)
def c_containing_cube(  np.ndarray[np.float64_t, ndim=1] atom_coordinate,
                        np.ndarray[np.float64_t, ndim=1] origin_crd,
                        np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                        np.ndarray[np.float64_t, ndim=1] spacing,
                        np.ndarray[np.int64_t, ndim=2]   eight_corner_shifts,
                        np.ndarray[np.float64_t, ndim=1] grid_x,
                        np.ndarray[np.float64_t, ndim=1] grid_y,
                        np.ndarray[np.float64_t, ndim=1] grid_z ):

    cdef:
        np.ndarray[np.float64_t, ndim=1] tmp
        np.ndarray[np.float64_t, ndim=1] corner_crd

        np.ndarray[np.int64_t, ndim=1]   lower_corner
        np.ndarray[np.int64_t, ndim=1]   corner

        list eight_corners, distances
        int nearest_ind, furthest_ind

    if not c_is_in_grid(atom_coordinate, origin_crd, uper_most_corner_crd):
        return [], 0, 0
    
    tmp = atom_coordinate - origin_crd
    lower_corner = np.array(tmp / spacing, dtype=int)
    eight_corners = [lower_corner + shift for shift in eight_corner_shifts]

    distances = []
    for corner in eight_corners:
        corner_crd = c_get_corner_crd(corner, grid_x, grid_y, grid_z)
        distances.append(cdistance(corner_crd, atom_coordinate))

    nearest_ind   = distances.index(min(distances))
    furthest_ind  = distances.index(max(distances))
    return eight_corners, nearest_ind, furthest_ind

@cython.boundscheck(False)
def c_lower_corner_of_containing_cube(  np.ndarray[np.float64_t, ndim=1] atom_coordinate,
                                        np.ndarray[np.float64_t, ndim=1] origin_crd,
                                        np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                                        np.ndarray[np.float64_t, ndim=1] spacing):
    cdef:
        np.ndarray[np.float64_t, ndim=1] tmp
        np.ndarray[np.int64_t, ndim=1]   lower_corner

    if not c_is_in_grid(atom_coordinate, origin_crd, uper_most_corner_crd):
        return np.array([], dtype=int)

    tmp = atom_coordinate - origin_crd
    lower_corner = np.array(tmp / spacing, dtype=int)
    return lower_corner

@cython.boundscheck(False)
def c_corners_within_radius(np.ndarray[np.float64_t, ndim=1] atom_coordinate,
                            double radius,
                            np.ndarray[np.float64_t, ndim=1] origin_crd,
                            np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                            np.ndarray[np.int64_t, ndim=1]   uper_most_corner,
                            np.ndarray[np.float64_t, ndim=1] spacing,
                            np.ndarray[np.float64_t, ndim=1] grid_x,
                            np.ndarray[np.float64_t, ndim=1] grid_y,
                            np.ndarray[np.float64_t, ndim=1] grid_z,
                            np.ndarray[np.int64_t, ndim=1]   gird_counts ):

    cdef:
        list corners
        int count_i, count_j, count_k
        int i, j, k
        float r, R 

        np.ndarray[np.int64_t, ndim=1] lower_corner
        np.ndarray[np.int64_t, ndim=1] corner

        np.ndarray[np.float64_t, ndim=1] lower_corner_crd
        np.ndarray[np.float64_t, ndim=1] corner_crd
        np.ndarray[np.float64_t, ndim=1] tmp
        np.ndarray[np.float64_t, ndim=1] lower_bound
        np.ndarray[np.float64_t, ndim=1] uper_bound
        np.ndarray[np.float64_t, ndim=1] dx2, dy2, dz2
        
    assert radius >= 0, "radius must be non-negative"
    if radius == 0:
        return []

    lower_corner = c_lower_corner_of_containing_cube(atom_coordinate, origin_crd, uper_most_corner_crd, spacing)
    if lower_corner.shape[0] > 0:

        lower_corner_crd = c_get_corner_crd(lower_corner, grid_x, grid_y, grid_z)
        r = radius + cdistance(lower_corner_crd, atom_coordinate)

        tmp = np.ceil(r / spacing)
        count_i, count_j, count_k = np.array(tmp, dtype=int)

        corners = []
        for i in range( -count_i, count_i + 1 ):
            for j in range( -count_j, count_j + 1 ):
                for k in range( -count_k, count_k + 1 ):

                    corner = lower_corner + np.array([i,j,k], dtype=int)

                    if np.all(corner >= 0) and np.all(corner <= uper_most_corner):
                        corner_crd = c_get_corner_crd(corner, grid_x, grid_y, grid_z)

                        if cdistance(corner_crd, atom_coordinate) <= radius:
                            corners.append(corner)
        return corners
    else:
        lower_bound = origin_crd - radius
        uper_bound  = uper_most_corner_crd + radius
        if np.any(atom_coordinate < lower_bound) or np.any(atom_coordinate > uper_bound):
            return []
        else:
            dx2 = (grid_x - atom_coordinate[0])**2
            dy2 = (grid_y - atom_coordinate[1])**2
            dz2 = (grid_z - atom_coordinate[2])**2

            corners = []
            count_i, count_j, count_k = gird_counts

            for i in range(count_i):
                for j in range(count_j):
                    for k in range(count_k):
                        R = dx2[i] + dy2[j] + dz2[k]
                        R = sqrt(R)
                        if R <= radius:
                            corners.append(np.array([i,j,k], dtype=int))
            return corners

@cython.boundscheck(False)
def c_is_row_in_matrix( np.ndarray[np.int64_t, ndim=1] row, 
                        list matrix):
    cdef:
        np.ndarray[np.int64_t, ndim=1] r

    for r in matrix:
        if (row == r).all():
            return True
    return False

@cython.boundscheck(False)
def c_ten_corners(  np.ndarray[np.float64_t, ndim=1] atom_coordinate,
                    np.ndarray[np.float64_t, ndim=1] origin_crd,
                    np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                    np.ndarray[np.int64_t, ndim=1]   uper_most_corner,
                    np.ndarray[np.float64_t, ndim=1] spacing,
                    np.ndarray[np.int64_t, ndim=2]   eight_corner_shifts,
                    np.ndarray[np.int64_t, ndim=2]   six_corner_shifts,
                    np.ndarray[np.float64_t, ndim=1] grid_x,
                    np.ndarray[np.float64_t, ndim=1] grid_y,
                    np.ndarray[np.float64_t, ndim=1] grid_z ):
    """
    to find the ten corners as described in the Qin et al J Chem Theory Comput 2014, 10, 2824
    """
    cdef:
        list eight_corners, six_corners, three_corners, ten_corners
        int nearest_ind, furthest_ind
        int i
        np.ndarray[np.int64_t, ndim=1] nearest_corner
        np.ndarray[np.int64_t, ndim=1] corner

    eight_corners, nearest_ind, furthest_ind = c_containing_cube(atom_coordinate, origin_crd,
                                                                uper_most_corner_crd, spacing,
                                                                eight_corner_shifts,
                                                                grid_x, grid_y, grid_z)
    if not eight_corners:
        raise RuntimeError("Atom is outside the grid")

    nearest_corner  = eight_corners[nearest_ind]
    for i in range(len(nearest_corner)):
        if nearest_corner[i] == 0 or nearest_corner[i] == uper_most_corner[i]:
            raise RuntimeError("The nearest corner is on the grid boundary")

    six_corners = [nearest_corner + corner for corner in six_corner_shifts]

    three_corners = []
    for corner in six_corners:
        if not c_is_row_in_matrix(corner, eight_corners):
            three_corners.append(corner)

    eight_corners.pop(furthest_ind)
    ten_corners = eight_corners + three_corners
    return ten_corners
    

@cython.boundscheck(False)
def c_distr_charge_one_atom( np.ndarray[np.float64_t, ndim=1] atom_coordinate, 
                             double charge,
                             np.ndarray[np.float64_t, ndim=1] origin_crd,
                             np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                             np.ndarray[np.int64_t, ndim=1]   uper_most_corner,
                             np.ndarray[np.float64_t, ndim=1] spacing,
                             np.ndarray[np.int64_t, ndim=2]   eight_corner_shifts,
                             np.ndarray[np.int64_t, ndim=2]   six_corner_shifts,
                             np.ndarray[np.float64_t, ndim=1] grid_x,
                             np.ndarray[np.float64_t, ndim=1] grid_y,
                             np.ndarray[np.float64_t, ndim=1] grid_z ):
    cdef:
        int i, j, k, row 
        list delta_vectors, ten_corners
        np.ndarray[np.int64_t, ndim=1] corner
        np.ndarray[np.float64_t, ndim=1] b_vector = np.zeros([10], dtype=float)
        np.ndarray[np.float64_t, ndim=2] a_matrix = np.zeros([10,10], dtype=float)
        np.ndarray[np.float64_t, ndim=1] corner_crd
        np.ndarray[np.float64_t, ndim=1] distributed_charges

    ten_corners = c_ten_corners(atom_coordinate, origin_crd, uper_most_corner_crd, uper_most_corner,
                                spacing, eight_corner_shifts, six_corner_shifts, grid_x, grid_y, grid_z)
    b_vector[0] = charge
    a_matrix[0,:] = 1.0

    delta_vectors = []
    for corner in ten_corners:
        corner_crd = c_get_corner_crd(corner, grid_x, grid_y, grid_z)
        delta_vectors.append(corner_crd - atom_coordinate)

    for j in range(10):
        a_matrix[1][j] = delta_vectors[j][0]
        a_matrix[2][j] = delta_vectors[j][1]
        a_matrix[3][j] = delta_vectors[j][2]

    row = 3
    for i in range(3):
        for j in range(i, 3):
            row += 1
            for k in range(10):
                a_matrix[row][k] = delta_vectors[k][i] * delta_vectors[k][j]

    distributed_charges = np.linalg.solve(a_matrix, b_vector)
    return ten_corners, distributed_charges


@cython.boundscheck(False)
def c_cal_potential_grid(   str name,
                            np.ndarray[np.float64_t, ndim=2] crd,
                            np.ndarray[np.float64_t, ndim=1] grid_x,
                            np.ndarray[np.float64_t, ndim=1] grid_y,
                            np.ndarray[np.float64_t, ndim=1] grid_z,
                            np.ndarray[np.float64_t, ndim=1] origin_crd,
                            np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                            np.ndarray[np.int64_t, ndim=1]   uper_most_corner,
                            np.ndarray[np.float64_t, ndim=1] spacing,
                            np.ndarray[np.int64_t, ndim=1]   gird_counts,
                            np.ndarray[np.float64_t, ndim=1] charges,
                            np.ndarray[np.float64_t, ndim=1] lj_sigma):

    cdef:
        list corners
        int natoms = crd.shape[0]
        int i_max = grid_x.shape[0]
        int j_max = grid_y.shape[0]
        int k_max = grid_z.shape[0]
        int i, j, k
        int atom_ind
        double charge, lj_diameter
        double d, exponent
        double dx_tmp, dy_tmp
        np.ndarray[np.float64_t, ndim=3] grid = np.zeros([i_max, j_max, k_max], dtype=np.float)
        np.ndarray[np.float64_t, ndim=3] grid_tmp
        np.ndarray[np.float64_t, ndim=1] atom_coordinate
        np.ndarray[np.float64_t, ndim=1] dx2, dy2, dz2

    if name != "occupancy":

        if name == "LJa":
            exponent = 3.
        elif name == "LJr":
            exponent = 6.
        elif name == "electrostatic":
            exponent = 0.5
        else:
            raise RuntimeError("Wrong grid name %s"%name)

        grid_tmp = np.empty([i_max, j_max, k_max], dtype=np.float)
        for atom_ind in range(natoms):
            atom_coordinate = crd[atom_ind]
            charge = charges[atom_ind]
            lj_diameter = lj_sigma[atom_ind]

            dx2 = (atom_coordinate[0] - grid_x)**2
            dy2 = (atom_coordinate[1] - grid_y)**2
            dz2 = (atom_coordinate[2] - grid_z)**2

            for i in range(i_max):
                dx_tmp = dx2[i]
                for j in range(j_max):
                    dy_tmp = dy2[j]
                    for k in range(k_max):

                        d = dx_tmp + dy_tmp + dz2[k]
                        d = d**exponent
                        grid_tmp[i,j,k] = charge / d

            corners = c_corners_within_radius(atom_coordinate, lj_diameter, origin_crd, uper_most_corner_crd,
                                                uper_most_corner, spacing, grid_x, grid_y, grid_z, gird_counts)

            for i, j, k in corners:
                grid_tmp[i,j,k] = 0.

            grid += grid_tmp
    else:
        for atom_ind in range(natoms):
            atom_coordinate = crd[atom_ind]
            lj_diameter = lj_sigma[atom_ind]
            corners = c_corners_within_radius(atom_coordinate, lj_diameter, origin_crd, uper_most_corner_crd,
                                                  uper_most_corner, spacing, grid_x, grid_y, grid_z, gird_counts)
            for i, j, k in corners:
                grid[i,j,k] = 1.
    return grid

@cython.boundscheck(False)
def c_cal_charge_grid(  str name,
                        np.ndarray[np.float64_t, ndim=2] crd,
                        np.ndarray[np.float64_t, ndim=1] charges,
                        np.ndarray[np.float64_t, ndim=1] origin_crd,
                        np.ndarray[np.float64_t, ndim=1] uper_most_corner_crd,
                        np.ndarray[np.int64_t, ndim=1]   uper_most_corner,
                        np.ndarray[np.float64_t, ndim=1] spacing,
                        np.ndarray[np.int64_t, ndim=2]   eight_corner_shifts,
                        np.ndarray[np.int64_t, ndim=2]   six_corner_shifts,
                        np.ndarray[np.float64_t, ndim=1] grid_x,
                        np.ndarray[np.float64_t, ndim=1] grid_y,
                        np.ndarray[np.float64_t, ndim=1] grid_z ):

    cdef:
        int atom_ind, i, l, m, n 
        int natoms = crd.shape[0]
        int i_max = grid_x.shape[0]
        int j_max = grid_y.shape[0]
        int k_max = grid_z.shape[0]
        double charge
        list ten_corners
        np.ndarray[np.float64_t, ndim=1] distributed_charges
        np.ndarray[np.float64_t, ndim=1] atom_coordinate
        np.ndarray[np.float64_t, ndim=3] grid = np.zeros([i_max, j_max, k_max], dtype=np.float)

    assert name in ["occupancy", "LJa", "LJr", "electrostatic"], "Name %s not allowed"%name

    if name != "occupancy":
        for atom_ind in range(natoms):
            atom_coordinate = crd[atom_ind]
            charge = charges[atom_ind]
            ten_corners, distributed_charges = c_distr_charge_one_atom( atom_coordinate, charge,
                                                                    origin_crd, uper_most_corner_crd,
                                                                    uper_most_corner, spacing,
                                                                    eight_corner_shifts, six_corner_shifts,
                                                                    grid_x, grid_y, grid_z)
            for i in range(len(ten_corners)):
                l, m, n = ten_corners[i]
                grid[l, m, n] += distributed_charges[i]
    else:
        for atom_ind in range(natoms):
            atom_coordinate = crd[atom_ind]
            ten_corners = c_ten_corners(atom_coordinate, origin_crd, uper_most_corner_crd,
                                    uper_most_corner, spacing, eight_corner_shifts, six_corner_shifts,
                                    grid_x, grid_y, grid_z )
            for i in range(len(ten_corners)):
                l, m, n = ten_corners[i]
                grid[l, m, n] = 1.0
    return grid


