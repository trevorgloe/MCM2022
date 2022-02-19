%% Testing distance calculation
close all;
clear;
clc;

%% distance()
lat = [45, 40, 39];
lon = [72, 72, 72];
dist = stdist(lat,lon);
disp(2*pi*6378*(dist/360)*2)

% [d1km d2km]=lldistkm(latlon1,latlon2)

latlon1 = [45, 72];
latlon2 = [40, 72];
latlon3 = [39, 72];

[d1kmA d2km] = lldistkm(latlon1,latlon2);
[d1kmB d2km] = lldistkm(latlon2,latlon3);
disp(d1kmA+d1kmB);

load('track_points2.csv');
T = readtable('track_points2.csv');
latlon = table2array(T);
dist2 = zeros(1,length(latlon)-1);
distcount = 0;
totaldist = zeros(1,length(latlon));

for i=1:length(latlon)-1
    dist2(i)=lldistkm([latlon(i,2), latlon(i,1)],[latlon(i+1,2), latlon(i+1,1)]);
    distcount = distcount+dist2(i);
    totaldist(i) = distcount;
end
totaldist(end) = totaldist(end-1);

disp("Loop answer: ")
disp(sum(dist2));

figure()
plot(totaldist,latlon(:,3));

dist3 = stdist(latlon(:,1),latlon(:,2));
disp(2*pi*6378*(dist3/360)*2)


% --Example 1, short distance:
%   latlon1=[-43 172];
%   latlon2=[-44  171];
%   [d1km d2km]=distance(latlon1,latlon2)
%   d1km =
%           137.365669065197 (km)
%   d2km =
%           137.368179013869 (km)
%   %d1km approximately equal to d2km

%arclen = distance('gc',[37,-76],[37,-9])