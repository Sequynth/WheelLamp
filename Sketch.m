% define two circles in space, rotated to each other by 90 degrees and
% their center offset from the origin by r/2 in each direction.

% number of points
N = 36;

% width of the felge
felgeWidth = 0.025; % m

% inner radius
r = 0.609/2; % (roughly a 28" wheel)

% get N equally spaced points
phi = linspace(0, 2*pi, N+1);
% we remove the last point, since it is equal to the first, also add a 90° rotation.
phi = phi(1:end-1) - pi/2;

% distance between both wheel centers
off = r * 2/3;

% position the points on a circle in space
w1(1 ,:) = r*sin(phi);
w1(2, :) = r*cos(phi);
w1(3, :) = 0;

% rotate a second set of points with respect to the ground
ang = 45;
w2 = (w1'*rotx(ang))';

% move the second set of points
w2(1, :) = w2(1, :) - off;

% rotate the numeration of the points by 180°
w2 = circshift(w2, N/2, 2);

% min height needed
height = r*sind(ang) + 0.1;

% lift both rings from the ground plate
w1(3, :) = w1(3, :) + height;
w2(3, :) = w2(3, :) + height;

%%

figure
hold on
plot3(w1(1,:),w1(2,:),w1(3,:),'ro')
plot3(w2(1,:),w2(2,:),w2(3,:),'bo')
plot3(w2(1,11),w2(2,11),w2(3,11),'x')
axis equal

% for ii = 1:N
%     plot3([w1(1,ii) w2(1,ii)],[w1(2,ii) w2(2,ii)],[w1(3,ii) w2(3,ii)])
% end

connections = zeros(N,2);

connections(:,1) = 1:N;
connections(2:2:N,2) = (2:2:N)';
connections(1:2:N,2) = [1; (N-1:-2:3)'];


connections(ceil(N/3)+2:floor(N*2/3),:) = [];

c = [25 13; 12 12; 13 25; 26 26; 27 11; 10 10; 11 27; 28 28;29 9; 8 8; 9 29;30 30; 31 7; 6 6; 7 31; 32 32; 33 5; 4 4; 5 33; 34 34; 35 3; 2 2; 3 35; 36 36; 1 1];

% for ii = 1:size(connections,1)
%     %plot3([w1(1,connections(ii,1)) w2(1,connections(ii,2))],[w1(2,connections(ii,1)) w2(2,connections(ii,2))],[w1(3,connections(ii,1)) w2(3,connections(ii,2))],'Color',hsv2rgb([ii/size(connections,1)/2 1 1]))
%     plot3([w1(1,c(ii,1)) w2(1,c(ii,2))],[w1(2,c(ii,1)) w2(2,c(ii,2))],[w1(3,c(ii,1)) w2(3,c(ii,2))],'Color',hsv2rgb([ii/size(c,1)/2 1 1]))
%     length(ii) = sqrt((w1(1,connections(ii,1))-w2(1,connections(ii,2))).^2 + (w1(2,connections(ii,1))-w2(2,connections(ii,2))).^2 + (w1(3,connections(ii,1))-w2(3,connections(ii,2))).^2);
% end

for ii = 1:4:21
    %plot3([w1(1,connections(ii,1)) w2(1,connections(ii,2))],[w1(2,connections(ii,1)) w2(2,connections(ii,2))],[w1(3,connections(ii,1)) w2(3,connections(ii,2))],'Color',hsv2rgb([ii/size(connections,1)/2 1 1]))
    plot3([w1(1,c(ii+1,1)) w2(1,c(ii+1,2))],[w1(2,c(ii+1,1)) w2(2,c(ii+1,2))],[w1(3,c(ii+1,1)) w2(3,c(ii+1,2))],'Color',hsv2rgb([ii/size(c,1)/2 1 1]))
    plot3([w1(1,c(ii+2,1)) w2(1,c(ii+2,2))],[w1(2,c(ii+2,1)) w2(2,c(ii+2,2))],[w1(3,c(ii+2,1)) w2(3,c(ii+2,2))],'Color',hsv2rgb([ii/size(c,1)/2 1 1]))
    plot3([w1(1,c(ii+3,1)) w2(1,c(ii+3,2))],[w1(2,c(ii+3,1)) w2(2,c(ii+3,2))],[w1(3,c(ii+3,1)) w2(3,c(ii+3,2))],'Color',hsv2rgb([ii/size(c,1)/2 1 1]))
    plot3([w1(1,c(ii+4,1)) w2(1,c(ii+4,2))],[w1(2,c(ii+4,1)) w2(2,c(ii+4,2))],[w1(3,c(ii+4,1)) w2(3,c(ii+4,2))],'Color',hsv2rgb([ii/size(c,1)/2 1 1]))
    length(ii) = sqrt((w1(1,connections(ii,1))-w2(1,connections(ii,2))).^2 + (w1(2,connections(ii,1))-w2(2,connections(ii,2))).^2 + (w1(3,connections(ii,1))-w2(3,connections(ii,2))).^2);
end



view([-23 24])

disp(connections)

%length of connections
fprintf('At least %.2f m of wire are necessary\n',sum(length));
fprintf('Approximately %.2f m are needed if one continuous wire is used\n',sum(length)+ 2*pi*r);

%% draw the shape of the wooden plate
margin = 0.05; % m
bottom  = -r - margin;
top     =  r + margin;
right   =  r + margin;
left    = -r - margin - off;

fprintf('necessary dimenions of the groundplate: %.1fmm x %.1fmm\n', (top - bottom)*1e3, (right - left)*1e3);

plot3([left right right left], [bottom bottom top top], [0 0 0 0]);

%% necessary posts

% red ring
% 4 times
Post1 = height + felgeWidth;

% blue Ring
% 1 time
Post2 = height -r*sind(ang)+2*felgeWidth;
% 1 times
Post3 = height +r*sind(ang)+2*felgeWidth;
% 2 times
Post4 = height +r*sind(ang)+2*felgeWidth;

Posts = [4 Post1;
    1 Post2;
    1 Post3;
    2 Post4]






