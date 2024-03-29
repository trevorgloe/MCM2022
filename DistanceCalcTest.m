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
inclination = zeros(1,length(latlon));

changeinelevation = zeros(1,length(latlon));

kmtft = 3280.84;

for i=1:length(latlon)-1
    dist2(i)=lldistkm([latlon(i,2), latlon(i,1)],[latlon(i+1,2), latlon(i+1,1)]);
    
    changeinelevation(i) = latlon(i+1,3)-latlon(i,3);
    inclination(i) = atand((latlon(i+1,3)-latlon(i,3))/(dist2(i)*kmtft));
    
    distcount = distcount+dist2(i);
    totaldist(i) = distcount;
end

iterations = round((length(latlon)-1)/5);
dist4 = zeros(1,iterations);
inclination2 = zeros(1,iterations);

distcount2 = 0;
totaldist2 = zeros(1,iterations);

for i=1:iterations-2
    dist4(i) = dist2(i*5)+dist2(i*5+1)+dist2(i*5+2)+dist2(i*5+3)+dist2(i*5+4);
    inclination2(i) = atand((latlon((i+1)*5,3)-latlon(i*5,3))/(dist4(i)*kmtft));
    
    distcount2 = distcount2+dist4(i);
    totaldist2(i) = distcount2;
end

inclination2 = inclination2*3;

totaldist(end) = totaldist(end-1);

disp("Max change in elevation [deg]: ");
disp(max(changeinelevation));
disp("Min change in elevation [deg]: ");
disp(min(changeinelevation));

disp("Max change in distance [ft]: ");
disp(max(dist2)*kmtft);
disp("Min change in distance [ft]: ");
disp(min(dist2)*kmtft);

disp("Max inclination1 [deg]: ");
disp(max(inclination));
disp("Min inclination1 [deg]: ");
disp(min(inclination));

disp("Max inclination2 [deg]: ");
disp(max(inclination2));
disp("Min inclination2 [deg]: ");
disp(min(inclination2));

inclination(end) = inclination(end-1);

inclinationgradientred = (abs(inclination/(abs(max(inclination)-min(inclination)))))*256;
inclinationgradientgreen = 256-(abs(inclination/(abs(max(inclination)-min(inclination)))))*256;
hsvcolor = zeros(3,length(inclinationgradientred))';
hsvcolortest = [2 3 4; 1 2 3; 4 5 6];
for i=1:length(inclinationgradientred)
    hsvcolor(i,1) = inclinationgradientred(i);
    hsvcolor(i,2) = 0;
    hsvcolor(i,3) = 0;
end

% hsvcolor = [inclinationgradientred; inclinationgradientgreen; zeros(1,length(inclinationgradientred))]';

disp("Loop answer: ")
disp(sum(dist2));

c = gradient(latlon(:,3));
% c = inclination;
figure()
scatter(totaldist,latlon(:,3),20,c);

figure()
plot(totaldist,inclination);

figure()
plot(totaldist2,inclination2*1.5);

figure()
scatter(latlon(:,1),latlon(:,2),10,inclination*4);
colorbar;

figure()
geoshow(latlon(:,2),latlon(:,1));
% geoshow(latlon(:,2),latlon(:,1),'Color',hsvcolor);

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