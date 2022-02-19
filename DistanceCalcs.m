%% Testing distance calculation
close all;
clear;
clc;

%% basic distance and inclination

load('flanders_track_points.csv');
T = readtable('flanders_track_points.csv');
latlon = table2array(T); % NOTE: column 3 in latlon is the elevation in ft

dist = zeros(1,length(latlon)-1); % distance between two latlon points
totaldist = zeros(1,length(latlon)); % total distance up to that point in the track
inclination = zeros(1,length(latlon));

kmtft = 3280.84; % km to ft conversion

for i=1:length(latlon)-1
    dist(i)=lldistkm([latlon(i,2), latlon(i,1)],[latlon(i+1,2), latlon(i+1,1)]); %lldistkm gives the distance between two latlon points in km, created by some guy on the internet
    
    % changeinelevation(i) = latlon(i+1,3)-latlon(i,3);
    inclination(i) = atand((latlon(i+1,3)-latlon(i,3))/(dist(i)*kmtft)); % inclination of the track given by change in elevation / distance (rise/run) in degrees
    % inclination calculated for single steps definitely undershoots the true grade provided by ridewithgps, so maybe multiply by 4?
    
    totaldist(i+1) = totaldist(i)+dist(i);
end

figure()
subplot(2,1,1)
plot(totaldist,inclination*4); % times 4 term is because shits wack bro
subplot(2,1,2)
plot(totaldist,latlon(:,3));

%% Calculation of smoother inclination for accuracy and nice graphs

iterations = round((length(latlon)-1)/5); % iterate for 1 per 5 data points
smooth_dist = zeros(1,iterations);
smooth_inclination = zeros(1,iterations);

smooth_totaldist = zeros(1,iterations);

for i=1:iterations-2
    smooth_dist(i) = dist(i*5)+dist(i*5+1)+dist(i*5+2)+dist(i*5+3)+dist(i*5+4);
    smooth_inclination(i) = atand((latlon((i+1)*5,3)-latlon(i*5,3))/(smooth_dist(i)*kmtft));
    
    smooth_totaldist(i+1) = smooth_totaldist(i)+smooth_dist(i);
end
smooth_totaldist(end) = smooth_totaldist(end-1);

figure()
subplot(3,1,1)
plot(totaldist,inclination*4); % times 4 term is because shits wack bro
subplot(3,1,2)
plot(smooth_totaldist,smooth_inclination*4); % times 4 term is because shits wack bro
subplot(3,1,3)
plot(totaldist,latlon(:,3));

disp("Total Distance (compare to ridewithgps) [mi]: ")
disp(sum(dist)*(kmtft/5280));

%% Plot a cool map of the track
c = gradient(latlon(:,3)); % color gradient that scatter() is able to understand for some reason, god bless

figure()
scatter(latlon(:,1),latlon(:,2),10,inclination*4);
colorbar;