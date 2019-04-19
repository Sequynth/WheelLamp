
import numpy as np
import matplotlib.pyplot as plt

# define two circles in space, rotated to each other by 90 degrees and
# their center offset from the origin by r/2 in each direction.

# number of points
N = 36;

# width of the felge
felgeWidth = 0.025; # m

# inner radius
r = 0.609/2; # (roughly a 28" wheel)

# get N equally spaced points
phi = linspace(0, 2*pi, N+1);
# we remove the last point, since it is equal to the first, also add a 90�
# rotation and half of the angle between two points.
phi = phi(1:end-1) - pi/2 + phi(2)/2;

# distance between both wheel centers
off = r * 2/3;

# position the points on a circle in space
w1(1 ,:) = r*sin(phi);
w1(2, :) = r*cos(phi);
w1(3, :) = 0;

# rotate a second set of points with respect to the ground
ang = 45;
w2 = (w1'*rotx(ang))';

# move the second set of points
w2(1, :) = w2(1, :) - off;

# rotate the numeration of the points by 180�
w2 = circshift(w2, N/2, 2);

# min height needed
height = r*sind(ang) + 0.1;

# lift both rings from the ground plate
w1(3, :) = w1(3, :) + height;
w2(3, :) = w2(3, :) + height;

# calculate the necessary support height at specified positions

## Show the wheels
# use matlab standard colors
colors = lines;

figure
hold on
for ii = 1:size(w1, 2)
    # annotate the points
    text(w1(1, ii), w1(2, ii), w1(3, ii), num2str(ii))
    text(w2(1, ii), w2(2, ii), w2(3, ii), num2str(ii))
end
axis equal

# set the connected pairs
#c = [25 13; 12 12; 13 25; 26 26; 27 11; 10 10; 11 27; 28 28;29 9; 8 8; 9 29;30 30; 31 7; 6 6; 7 31; 32 32; 33 5; 4 4; 5 33; 34 34; 35 3; 2 2; 3 35; 36 36; 1 1];

c = [25 12; 11 11; 12 25; 26 26;...
    27 10; 9 9; 10 27; 28 28;...
    29 8; 7 7; 8 29; 30 30;...
    31 6; 5 5; 6 31; 32 32;...
    33 4; 3 3; 4 33; 34 34;...
    35 2; 1 1; 2 35; 36 36];
# numbers on wheel one
c1 = c(:, 1);
# and two
c2 = c(:, 2);

# how many lines from one single rope
n = 4;

idx = 1;
for ii = 1:4:21
    xcoords = [w1(1,c1(ii)) w2(1,c2(ii)) w2(1,c2(ii+1)) w1(1,c1(ii+1)) w1(1,c1(ii+2)) w2(1,c2(ii+2)) w2(1,c2(ii+3)) w1(1,c1(ii+3)) w1(1,c1(ii))];
    ycoords = [w1(2,c1(ii)) w2(2,c2(ii)) w2(2,c2(ii+1)) w1(2,c1(ii+1)) w1(2,c1(ii+2)) w2(2,c2(ii+2)) w2(2,c2(ii+3)) w1(2,c1(ii+3)) w1(2,c1(ii))];
    zcoords = [w1(3,c1(ii)) w2(3,c2(ii)) w2(3,c2(ii+1)) w1(3,c1(ii+1)) w1(3,c1(ii+2)) w2(3,c2(ii+2)) w2(3,c2(ii+3)) w1(3,c1(ii+3)) w1(3,c1(ii))];
    plot3(xcoords, ycoords, zcoords, 'Color', colors(idx, :));    
    lengths(idx) = cumDist(xcoords, ycoords, zcoords);
    idx = idx+1;
end

view([-23 24])


#length of connections
fprintf('At least #.2f m of wire are necessary\n',sum(lengths));
lengths'

## draw the shape of the wooden plate
margin = 0.05; # m
bottom  = -r - margin;
top     =  r + margin;
right   =  r + margin;
left    = -r - margin - off;

fprintf('necessary dimenions of the groundplate: #.1fmm x #.1fmm\n', (top - bottom)*1e3, (right - left)*1e3);

plot3([left right right left], [bottom bottom top top], [0 0 0 0]);

## necessary posts

# red ring
# 4 times
Post1 = height + felgeWidth;

# blue Ring
# 1 time
Post2 = height -r*sind(ang)+2*felgeWidth;
# 1 times
Post3 = height +r*sind(ang)+2*felgeWidth;
# 2 times
Post4 = height +r*sind(ang)+2*felgeWidth;

Posts = [4 Post1;
    1 Post2;
    1 Post3;
    2 Post4]






