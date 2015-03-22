import sys
import numpy as np


def xrot_matrix(anglex, degree=True):
    if degree:
        anglex = np.pi * anglex / 180.0
    Rx = np.array([[1,       0,               0       ],
                   [0, np.cos(anglex), -np.sin(anglex)],
                   [0, np.sin(anglex),  np.cos(anglex)]])
    return Rx

def yrot_matrix(angley, degree=True):
    if degree:
        angley = np.pi * angley / 180.0
    Ry = np.array([[ np.cos(angley), 0, np.sin(angley)],
                   [       0,        1,       0       ],
                   [-np.sin(angley), 0, np.cos(angley)]]) 
    return Ry

def zrot_matrix(anglez, degree=True):
    if degree:
        anglez = np.pi * anglez / 180.0
    Rz = np.array([[np.cos(anglez), -np.sin(anglez), 0],
                   [np.sin(anglez),  np.cos(anglez), 0],
                   [      0,               0,        1]])
    return Rz
 

def rotation_matrix(axes, angles, degree=True):
    """return the 3D rotation matrix around the given angles"""

    R = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [ 0, 0, 1]])

    for axis, angle in zip(axes, angles):
        if axis == "x":
            R = np.dot(xrot_matrix(angle, degree), R)
        elif axis == "y":
            R = np.dot(yrot_matrix(angle, degree), R)
        elif axis == "z":
            R = np.dot(zrot_matrix(angle, degree), R)
        else:
            print "rotation axis can only be x, y or z."
            sys.exit(1)

    return R


def rotate(vector, axes, angles, degree=True):
    R = rotation_matrix(axes, angles, degree)
    return np.dot(R, vector.reshape(3))


def rotate_multi(vector, axes, angles, degree=True):
    R = rotation_matrix(axes, angles, degree)
    result = np.zeros_like(vector)
    for i in range(vector.size/3):
        result[i*3:i*3+3] = np.dot(R, vector[i*3:i*3+3])
    return result

def shift(vector, offsetvector):
    return vector + offsetvector


def shift_random(vector):
    return vector + np.random.rand(vector.shape)
