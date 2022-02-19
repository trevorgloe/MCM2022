load(course.data)
alt = table2array(T(:,3));
dist = dist2;

% dalt(ii) = 0;
for ii = 2:length(alt)-1
    dalt(ii-1) = alt(ii) - alt(ii-1);
end

n=50;
ii=1;
for jj = 1:length(dalt)
    if mod(jj,n) == 0
        dalt_s(ii) = abs(dalt(jj));
        dist_s(ii) = dist(jj);
        ii = ii + 1;
    end
end

phi = atand(dalt_s./dist_s);
course.phi = phi';