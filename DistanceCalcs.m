%% Testing distance calculation
% close all;
% clear;
% clc;

%% basic distance and inclination

load('tokyo_track_points.csv');
T = readtable('tokyo_track_points.csv');
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

%% Calculate the radius of curvature
% tangentvec = zeros(2,length(latlon))';
% radius = zeros(1,length(tangentvec))';
% 
% for i=1:length(latlon)-1
%     tangentvec(i,1) = latlon(i+1,1)-latlon(i,1);
%     tangentvec(i,2) = latlon(i+1,2)-latlon(i,2);
%     
%     radius(i) = abs(1/(norm(tangentvec(i+1,:))-norm(tangentvec(i,:))));
% end
% radius(end) = radius(end-1);
% disp("Minimum Radius");
% disp(min(abs(radius)));
% 
% [~, ind] = min(abs(radius));
% disp("Index of Minimum Radius");
% disp(ind);
% 
% totaldistmi = totaldist*(kmtft/5280);
% 
% % radius_nonzero = nonzeros(radius);
% % disp("Minimum Radius");
% % disp(min(abs(radius_nonzero)));
% 
% tangent(end) = tangent(end-1);
% radius(end) = radius(end-1);
% tangentvec(end,1) = tangentvec(end-1,1); 
% tangentvec(end,2) = tangentvec(end-1,2);
% 
% tangentvec = tangentvec.*kmtft; % convert to ft
% 
% figure()
% subplot(3,1,1)
% plot(totaldist,inclination*4); % times 4 term is because shits wack bro
% subplot(3,1,2)
% plot(totaldist,radius);
% subplot(3,1,3)
% plot(totaldist,latlon(:,3));
% sgtitle("Radius of curvature");

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

% disp("Total Distance (compare to ridewithgps) [mi]: ")
% disp(sum(dist)*(kmtft/5280));

disp("Total Distance (compare to ridewithgps) [mi]: ")
disp(totaldist(end)*(kmtft/5280));

% disp("Total Smooth Distance (compare to ridewithgps) [mi]: ")
% disp(sum(smooth_dist)*(kmtft/5280));

disp("Total Smooth Distance (compare to ridewithgps) [mi]: ")
disp(smooth_totaldist(end)*(kmtft/5280));

%% Plot a cool map of the track

c = gradient(latlon(:,3)); % color gradient that scatter() is able to understand for some reason, god bless

figure()
scatter(latlon(:,1),latlon(:,2),10,inclination*4);
colorbar;
title("Raw Data: ");

%% Sampled data track (221 points instead of 1105)
latlon_5sample = latlon(1:5:end,:); % NOTE: column 3 in latlon is the elevation in ft
inclination_5sample = inclination(:,1:5:end);
totaldist_5sample = totaldist(:,1:5:end);
% smooth_totaldist_5sample = smooth_totaldist(:,1:5:end);

tangent_5sample = zeros(1,length(latlon_5sample));
radius_5sample = zeros(1,length(tangent_5sample));

for i=1:length(latlon_5sample)-1
    tangent_5sample(i) = (latlon_5sample(i+1,2)-latlon_5sample(i,2))/(latlon_5sample(i+1,1)-latlon_5sample(i,1));
    radius_5sample(i) = 1/tangent_5sample(i);
end
tangent_5sample(end) = tangent_5sample(end-1);
radius_5sample(end) = radius_5sample(end-1);

disp("Index of max radius");
[~, ind] = max(radius_5sample);
disp(ind);

totaldistmi_5sample = totaldist_5sample*(kmtft/5280); % total distance up to that point in the track

figure()
scatter(latlon_5sample(:,1),latlon_5sample(:,2),10,inclination_5sample*4);
colorbar;
title("Sampled Data: ");

figure()
subplot(3,1,1)
plot(totaldist_5sample,inclination_5sample*4); % times 4 term is because shits wack bro
% subplot(3,1,2)
% plot(smooth_totaldist_5sample,smooth_inclination*4); % times 4 term is because shits wack bro
subplot(3,1,2)
plot(totaldist_5sample,radius_5sample); % times 4 term is because shits wack bro
subplot(3,1,3)
plot(totaldist_5sample,latlon_5sample(:,3));
sgtitle("Sampled Data");

% %% Interpolation of the smooth inclination [DOESN'T WORK]
% 
% % unique_figureout = unique(smooth_totaldist);
% [~, ind] = unique(smooth_totaldist);
% F = griddedInterpolant(smooth_totaldist(ind),smooth_inclination(ind));
% 
% interxq = linspace(0,smooth_totaldist(end),length(totaldist));
% intervq = F(interxq);
% 
% figure()
% plot(smooth_totaldist,smooth_inclination,'ro')
% hold on
% plot(interxq,intervq,'.')
% legend('Sample Points','Interpolated Values')
% 
% disp("Mean inclination: ");
% disp(mean(vq));
% % inter_dist = zeros(1,length(latlon)-1); % distance between two latlon points
% % inter_totaldist = zeros(1,length(latlon)); % total distance up to that point in the track
% % inter_inclination = zeros(1,length(latlon));
% % 
% % % inter_dist = 
% % for i=1:length(latlon)-5
% %     % inter_dist(i) = smooth_dist(ceil(i/5))+mod(i,5)*(smooth_dist(ceil(i/5)+1)-smooth_dist(ceil(i/5)));
% %     % inter_dist(i) = smooth_dist(ceil(i/5))+mod(i,5)*(smooth_dist(ceil(i/5)+1)-smooth_dist(ceil(i/5)));
% %     % inter_inclination(i) = atand((latlon(i+1,3)-latlon(i,3))/(inter_dist(i)*kmtft));
% %     
% %     inter_inclination(i) = smooth_inclination(ceil(i/5))+(mod(i,5)/5)*(smooth_inclination(ceil(i/5)+1)-smooth_inclination(ceil(i/5)));
% %     
% %     % inter_totaldist(i+1) = inter_totaldist(i)+inter_dist(i);
% %     
% %     inter_totaldist(i+1) = inter_totaldist(i)+(mod(i,5)/5)*(smooth_dist(ceil(i/5)+1)-smooth_dist(ceil(i/5)));
% %     
% %     % changeinelevation(i) = latlon(i+1,3)-latlon(i,3);
% %     % inclination(i) = atand((latlon(i+1,3)-latlon(i,3))/(dist(i)*kmtft)); % inclination of the track given by change in elevation / distance (rise/run) in degrees
% %     % inclination calculated for single steps definitely undershoots the true grade provided by ridewithgps, so maybe multiply by 4?
% % end
% % 
% % disp("Total Interpolated Distance (compare to ridewithgps) [mi]: ")
% % disp(sum(inter_dist)*(kmtft/5280));
% % 
% % disp("Total Interpolated Distance variable (compare to ridewithgps) [mi]: ")
% % disp(inter_totaldist(end)*(kmtft/5280));
% % 
% figure()
% subplot(4,1,1)
% plot(totaldist,inclination*4); % times 4 term is because shits wack bro
% subplot(4,1,2)
% plot(smooth_totaldist,smooth_inclination*4); % times 4 term is because shits wack bro
% subplot(4,1,3)
% plot(xq,vq*4); % times 4 term is because shits wack bro
% subplot(4,1,4)
% plot(totaldist,latlon(:,3));
% 
% %% Plotting the map with Power
% F = griddedInterpolant(x,P);
% 
% T = readtable('flanders_track_points.csv');
% latlon = table2array(T); % NOTE: column 3 in latlon is the elevation in ft
% 
% xq = linspace(0,x(end),length(latlon));
% vq = F(xq)';
% 
% % c = gradient(vq); % color gradient that scatter() is able to understand for some reason, god bless
% 
% 
% figure()
% subplot(1,2,1)
% scatter(latlon(:,1),latlon(:,2),10,intervq);
% colormap jet
% colorbar;
% title("Grade Data: ");
% 
% subplot(1,2,2)
% scatter(latlon(:,1),latlon(:,2),10,vq);
% colormap jet
% colorbar;
% title("Power Data: ");