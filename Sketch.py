import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math

def rotx(ang):
    # ang: rotation angle around the x-axis
    return np.array([[1, 0, 0], [0, math.cos(ang), -math.sin(ang)], [0, math.sin(ang), math.cos(ang)]])

def cumDist(x, y, z):
    # x: x-positions
    # y: y-positions
    # z: z-positions
    distance = 0;
    for ii in range(0, len(x)-1):
        distance = distance + math.sqrt( (x[ii]-x[ii+1])**2 + (y[ii]-y[ii+1])**2 + (z[ii]-z[ii+1])**2)
    return distance

def getCoords(w1, w2, c1, c2):
    # w1: coordinates on ring 1
    # w2: coordinates on ring 2
    # c1: holes on ring 1
    # c2: holes on ring 2
    pts = sum(~np.isnan(c1))
    coords = np.zeros(shape=(2*pts+1))
    for kk in range(0, pts, 1):
        if kk % 2 == 0:
            coords[2*kk]   = w1[int(c1[kk])]
            coords[2*kk+1] = w2[int(c2[kk])]
        else:
            coords[2*kk]   = w2[int(c2[kk])]
            coords[2*kk+1] = w1[int(c1[kk])]
    coords[2*pts] = w1[int(c1[0])]
    return coords

# define two circles in space, rotated to each other by 90 degrees and
# their center offset from the origin by r/2 in each direction.

# number of points
N = 36;

# inner radius
r = 0.611/2; # (roughly a 28" wheel)
# outer radius
r2 = 0.636/2; # m
# width of the felge
felgeWidth = r2-r; # m

# get N equally spaced points
phi = np.linspace(0, 2*math.pi, N+1);
# we remove the last point, since it is equal to the first, also add a 90�
# rotation and half of the angle between two points.
phi = phi[:-1] - math.pi/2 + phi[1]/2;

# distance between both wheel centers
off = r * 2/3;

# position the points on a circle in space
w1 = np.zeros(shape=(3, N))
w1[0, :] = r*np.sin(phi);
w1[1, :] = r*np.cos(phi);

# rotate a second set of points with respect to the ground
ang = 45;
rotmat = rotx(ang)
w2 = np.zeros(shape=(3, N))
for ii in range(1, N):
    w2[:, ii] = rotmat.dot(w1[:, ii])

# move the second set of points
w2[0, :] = w2[0, :] - off;
print('The rings are shifted by %.4f m' % off)

# rotate the numeration of the points by 180deg
w2 = np.roll(w2, int(N/2), 1);

# min height needed
height = r*math.sin(math.radians(ang)) + 0.1;
print('The rings are lifted by %.4f m' % height)
# lift both rings from the ground plate
w1[2, :] = w1[2, :] + height;
w2[2, :] = w2[2, :] + height;

# calculate the necessary support height at specified positions

## Show the wheels
# use matlab standard colors
# colors = lines;

fig = plt.figure()
ax = fig.gca(projection = '3d', proj_type = 'ortho')
for ii in range(0, w1.shape[1]):
    # annotate the points
    ax.text(w1[0, ii], w1[1, ii], w1[2, ii], ii)
    ax.text(w2[0, ii], w2[1, ii], w2[2, ii], ii)


# set the connected pairs
#c = [25 13; 12 12; 13 25; 26 26; 27 11; 10 10; 11 27; 28 28;29 9; 8 8; 9 29;30 30; 31 7; 6 6; 7 31; 32 32; 33 5; 4 4; 5 33; 34 34; 35 3; 2 2; 3 35; 36 36; 1 1];

c = np.array([[[24, 11], [10, 10], [11, 24], [25, 25]],
    [[26, 9], [8, 8], [9, 26], [27, 27]],
    [[28, 7], [6, 6], [7, 28], [29, 29]],
    [[30, 5], [4, 4], [5, 30], [31, 31]],
    [[32, 3], [2, 2], [3, 32], [33, 33]],
    [[34, 1], [0, 0], [1, 34], [35, 35]]])

# fist column: flat ring
# second column: rotated ring
c = np.array([[[24, 11], [10, 10], [11, 24], [25, 25]],
    [[26, 31], [27, 30], [np.nan, np.nan], [np.nan, np.nan]],
    [[ 8,  5], [ 9,  4], [np.nan, np.nan], [np.nan, np.nan]],
    [[ 2, 29], [ 3, 28], [np.nan, np.nan], [np.nan, np.nan]],
    [[32,  7], [33,  6], [np.nan, np.nan], [np.nan, np.nan]],
    [[30,  3], [31,  2], [np.nan, np.nan], [np.nan, np.nan]],
    [[ 5, 32], [ 4, 33], [np.nan, np.nan], [np.nan, np.nan]],
    [[34,  1], [ 0,  0], [ 1, 34], [35, 35]]])

c = np.array([[[22, 10], [24, 11], [np.nan, np.nan], [np.nan, np.nan]],
    [[11, 24], [13, 25], [np.nan, np.nan], [np.nan, np.nan]],
    [[25, 31], [26, 29], [np.nan, np.nan], [np.nan, np.nan]],
    [[ 9,  6], [10,  4], [np.nan, np.nan], [np.nan, np.nan]],
    [[ 2, 28], [ 3, 26], [np.nan, np.nan], [np.nan, np.nan]],
    [[32,  9], [33,  7], [np.nan, np.nan], [np.nan, np.nan]],
    [[30,  3], [31,  1], [np.nan, np.nan], [np.nan, np.nan]],
    [[ 5, 32], [ 4, 34], [np.nan, np.nan], [np.nan, np.nan]],
    [[34,  0], [ 1, 35], [np.nan, np.nan], [np.nan, np.nan]]])
# numbers on wheel one
c1 = c[:, :, 0];
# and two
c2 = c[:, :, 1];

# number of ropes
nRopes = c.shape[0]
# number of holes on each wheel
nHoles = c.shape[1]

coords = np.zeros(shape=(3, nRopes, 2*nHoles+1))
lengths = np.zeros(nRopes)
idx = 0;
for ii in range(0, nRopes, 1):
    # loop over ropes
    nPtsRope = sum(~np.isnan(c1[ii]))
    for jj in range(0, 3, 1):
        # loop over x, y and z
        coords[jj, ii, :(2*nPtsRope+1)] = getCoords(w1[jj, :], w2[jj, :], c1[ii, :], c2[ii, :])
    ax.plot(coords[0, ii, :(2*nPtsRope+1)], coords[1, ii, :(2*nPtsRope+1)], coords[2, ii, :(2*nPtsRope+1)]);
    lengths[idx] = cumDist(coords[0, ii, :(2*nPtsRope+1)], coords[1, ii, :(2*nPtsRope+1)], coords[2, ii, :(2*nPtsRope+1)]);
    idx = idx+1;

#ax.view_init(-23, 24)

#length of connections
print('At least %.2f m of wire are necessary' % (sum(lengths)));
print(lengths)

## draw the shape of the wooden plate
margin = 0.05; # m
bottom  = -r - margin;
top     =  r + margin;
right   =  r + margin;
left    = -r - margin - off;

print('necessary dimenions of the groundplate: %.1fmm x %.1fmm' % ((top - bottom)*1e3, (right - left)*1e3));

ax.plot([left, right, right, left], [bottom, bottom, top, top], [0, 0, 0, 0]);

## necessary posts

# red ring
# 4 times
Post1 = height + felgeWidth;

# blue Ring
# 1 time
Post2 = height -r*math.sin(math.radians(ang))+2*felgeWidth;
# 1 times
Post3 = height +r*math.sin(math.radians(ang))+2*felgeWidth;
# 2 times
Post4 = height +r*math.sin(math.radians(ang))+2*felgeWidth;

Posts = np.array([[4, Post1],
    [1, Post2],
    [1, Post3],
    [2, Post4]])

plt.show()
#plt.axis('equal')
