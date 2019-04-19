function distance = cumDist(x, y, z)
% calculate the cumulated distance between multiple points in 3d

% make sure input is of correct size
if ~isequal(size(x), size(y), size(z))
    error('input mus be of equal size')
end

% and a vector
if ~(isvector(x) && isvector(y) && isvector(z))
    error('input must be vectors')
end

distance = 0;
for ii = 1:numel(x)-1
    distance = distance + sqrt((x(ii)-x(ii+1))^2 + (y(ii)-y(ii+1))^2 + (z(ii)-z(ii+1))^2);
end