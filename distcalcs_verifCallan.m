%% Get distance data

T = readtable('flanders_track_points.csv');
latlon = table2array(T); % NOTE: column 3 in latlon is the elevation in ft
alt = 0.3048*table2array(T(:,3)); %[m]
lat = table2array(T(:,1));
long = table2array(T(:,2)); %deg

%% Lat long to Cartesian
[xEast,yNorth] = latlon2local(lat,long,alt,[lat(1),long(1),alt(1)]);

figure;
scatter(xEast,yNorth,10)
axis('equal'); % set 1:1 aspect ratio to see real-world shape

% R = 6378e3;
% x = cosd(lat).*cosd(long)*R;
% y = cosd(lat).*sind(long)*R;

[L,R,k] = curvature([xEast yNorth]);
format shortG
disp(L(end)/1000)
disp(xEast');

figure();
plot(L,R)
title('Curvature radius vs. cumulative curve length')
xlabel L
ylabel R
ylim([0,1000]);

figure()
h = plot(xEast,yNorth); grid on; axis equal
set(h,'marker','.');
xlabel x
ylabel y
title('2D curve with curvature vectors')
hold on
quiver(xEast,yNorth,k(:,1),k(:,2));
hold off

%% Functions
function [L,R,k] = curvature(X)
    % Radius of curvature and curvature vector for 2D or 3D curve
    %  [L,R,k] = curvature(X)
    %   X:   2 or 3 column array of x, y (and possibly z) coordiates
    %   L:   Cumulative arc length
    %   R:   Radius of curvature
    %   k:   Curvature vector
    % The scalar curvature value is 1./R
    % Version 2.6: Calculates end point values for closed curve
      N = size(X,1);
      dims = size(X,2);
      if dims == 2
        X = [X,zeros(N,1)];  % Use 3D expressions for 2D as well
      end
      L = zeros(N,1);
      R = NaN(N,1);
      k = NaN(N,3);
      for i = 2:N-1
        [R(i),~,k(i,:)] = circumcenter(X(i,:)',X(i-1,:)',X(i+1,:)');
        L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
      end
      if norm(X(1,:)-X(end,:)) < 1e-10 % Closed curve. 
        [R(1),~,k(1,:)] = circumcenter(X(end-1,:)',X(1,:)',X(2,:)');
        R(end) = R(1);
        k(end,:) = k(1,:);
        L(end) = L(end-1) + norm(X(end,:)-X(end-1,:));
      end
      i = N;
      L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
      if dims == 2
        k = k(:,1:2);
      end
end

function [R,M,k] = circumcenter(A,B,C)
% Center and radius of the circumscribed circle for the triangle ABC
%  A,B,C  3D coordinate vectors for the triangle corners
%  R      Radius
%  M      3D coordinate vector for the center
%  k      Vector of length 1/R in the direction from A towards M
%         (Curvature vector)
  D = cross(B-A,C-A);
  b = norm(A-C);
  c = norm(A-B);
  if nargout == 1
    a = norm(B-C);     % slightly faster if only R is required
    R = a*b*c/2/norm(D);
    if norm(D) == 0
      R = Inf;
    end
    return
  end
  E = cross(D,B-A);
  F = cross(D,C-A); 
  G = (b^2*E-c^2*F)/norm(D)^2/2;
  M = A + G;
  R = norm(G);  % Radius of curvature
  if R == 0
    k = G;
  elseif norm(D) == 0
    R = Inf;
    k = D;
  else
    k = G'/R^2;   % Curvature vector
  end
end

