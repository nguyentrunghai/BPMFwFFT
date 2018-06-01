
import numpy as np
import netCDF4

from IO import InpcrdLoad


def _rotation_matrix(u):
    """
    taken from AlGDock.Integrators.ExternalMC
    """
    u = np.array(u, dtype=float)
    assert np.all(0 <= u) and np.all( u<=1 ) , "u must be in between 0 and 1"
    # quaternion
    q = np.array( [ np.sqrt(1 - u[0]) * np.sin(2 * np.pi * u[1]),
                    np.sqrt(1 - u[0]) * np.cos(2 * np.pi * u[1]),
                    np.sqrt(u[0]) * np.sin(2 * np.pi * u[2]),
                    np.sqrt(u[0]) * np.cos(2 * np.pi * u[2]) ] )

    # Convert the quaternion into a rotation matrix
    rotMat = np.zeros([3,3], dtype=float)

    rotMat[0,0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]
    rotMat[0,1] = 2*q[1]*q[2] - 2*q[0]*q[3]
    rotMat[0,2] = 2*q[1]*q[3] + 2*q[0]*q[2]

    rotMat[1,0] = 2*q[1]*q[2] + 2*q[0]*q[3]
    rotMat[1,1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]
    rotMat[1,2] = 2*q[2]*q[3] - 2*q[0]*q[1]

    rotMat[2,0] = 2*q[1]*q[3] - 2*q[0]*q[2]
    rotMat[2,1] = 2*q[2]*q[3] + 2*q[0]*q[1]
    rotMat[2,2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]

    return rotMat


def _random_rotation_matrix():
    u = np.random.uniform(size=3)
    return _rotation_matrix(u)


def _rotate_molecule(A, crd):
    """
    A:  3X3 array, the rotation matrix
    crd:    natomsX3 array, coordinates of atoms in the molecule
    """
    assert A.shape == (3,3) and crd.shape[1] == 3, "A or crd's shape is wrong" 
    rotated_crd = np.empty(crd.shape, dtype=float)
    for i in range(crd.shape[0]):
        rotated_crd[i] = np.dot(A, crd[i])
    return rotated_crd


def _open_nc(crd, nc_file_name):
    nc_handle = netCDF4.Dataset(nc_file_name, "w", format="NETCDF4")
    nc_handle.createDimension("natoms", crd.shape[0])
    nc_handle.createDimension("Cartesian", crd.shape[1])
    nc_handle.createDimension("nrotations", None)
    nc_handle.createVariable("positions", "f4", ("nrotations", "natoms", "Cartesian"))
    return nc_handle


def _move_molecule_to_origin(crd):
    ligand_center = np.array([crd[:,i].mean() for i in range(3)], dtype=float)
    displacement  = -ligand_center
    for atom_ind in range(len(crd)):
        crd[atom_ind] += displacement
    return crd


def systematic_gen_rotation(ligand_crd, total_count, output_nc):
    """
    ligand_crd:    str, file name
    count:          int, total number of rotations
    output_nc:      str, name of output file
    """
    initial_crd = InpcrdLoad(ligand_crd).get_coordinates()
    initial_crd = _move_molecule_to_origin(initial_crd)

    total_count = np.ceil( total_count**(1./3) )
    u = np.linspace(0., 1., total_count+2)
    u = u[1:-1]
    print("Generating total %d rotations"%(len(u)**3))
    nc_handle = _open_nc(initial_crd, output_nc)

    count = -1
    for i in u:
        for j in u:
            for k in u:

                count += 1
                if count%100 == 0:
                    print("Doing %d-th"%count)

                A = _rotation_matrix([i, j, k])
                crd = _rotate_molecule(A, initial_crd)
                nc_handle.variables["positions"][count,:,:] = crd
    nc_handle.close()
    print("Generation Done!")
    return None


def random_gen_rotation(ligand_crd, total_count, output_nc):
    """
    ligand_crd:    str, file name
    total_count:    int, total number of rotations
    output_nc:      str, name of output file
    """
    initial_crd = InpcrdLoad(ligand_crd).get_coordinates()
    initial_crd = _move_molecule_to_origin(initial_crd)

    print("Generating total %d rotations"%total_count)
    nc_handle = _open_nc(initial_crd, output_nc)
    nc_handle.variables["positions"][0,:,:] = initial_crd

    for i in range(total_count):
        if (i%100) == 0:
            print("Doing %d-th"%i)
        A = _random_rotation_matrix()
        crd = _rotate_molecule(A, initial_crd)
        nc_handle.variables["positions"][i+1,:,:] = crd
    nc_handle.close()
    print("Generation Done!")
    return None


def random_rotation(inpcrd):
    """
    inpcrd: np.array of shape (natoms, 3)
    """
    initial_crd = _move_molecule_to_origin(inpcrd)
    A = _random_rotation_matrix()
    crd = _rotate_molecule(A, initial_crd)
    return crd

