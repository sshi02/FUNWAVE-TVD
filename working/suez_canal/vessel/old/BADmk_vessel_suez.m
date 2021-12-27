% PRODUCE VESSEL FILE FOR suez_canal

% - CREDIT -------------------------------------------------------------- 
% Maritime Casuality Specialist. (2021). Ever Given grounding in 
% Suez Canal - AIS based Dynamic Reconstruction. Retrieved July
% 7, 2021 at vimeo.com/531626438
%
% VesselFinder. (2021). EVER GIVEN, Container Ship - Details and current 
% position - IMO 9811000 MMSI 353136000 - VesselFinder. Retrieved July
% 7, 2021 at www.vesselfinder.com/vessels/EVER-GIVEN-IMO-9811000-MMSI-353136000

% - DESCRIPTION----------------------------------------------------------
% The following file will produce a vessel file of the Ever Given, with the
% following dimensions 400 x 59 m and a draught of 15 m. 
% The positional data will be generated from vessel_data.csv. Due to 
% complications between dense data and FUNWAVE, the positional data will 
% be smoothed out via regression. 

clear all;

vlf = fopen('vessel_00001', 'w');

fprintf(vlf, 'Title: Vessel # 1  - Ever Given\n');
fprintf(vlf, 'PRESSURE, 1\n');
fprintf(vlf, ['Length(m), Width(m), Alpha1(m),  ' ...
    'Alpha2(m), Beta(m), P(unit)\n']);
fprintf(vlf, '400.0  59.0, 0.5, 0.5, 0.5, 15.0\n');
fprintf(vlf, 'Time, X(m), Y(m)\n');

% GENERATING POSITION DATA
filename = 'vessel_data.csv';
vdata = csvread(filename, 1, 0);
x = 200; % BAD/UNVERIFIED INITIAL POS
y = 220; % 204, 267

% INTIAL COURSE
x1 = 1:1:36;
y1 = 1:1:36;
for ind = 1:36
    x1(ind) = x;
    y1(ind) = y;
    v = (vdata(ind, 2) + vdata(ind + 1, 2)) * 0.25;
    dt = vdata(ind + 1, 1) - vdata(ind, 1);
    x = x + dt * v * cosd(81 - (vdata(ind, 4) + vdata(ind + 1, 4)) * 0.5);
    y = y + dt * v * sind(81 - (vdata(ind, 4) + vdata(ind + 1, 4)) * 0.5);
end

x1 = x1(:);
y1 = y1(:);
X1 = [ones(size(x1)) x1];
b1 = regress(y1, X1);

% FIRST DIRECTIONAL CHANGE
x2 = 30:1:75;
y2 = 30:1:75;
for ind = 20:75
    x2(ind - 19) = x;
    y2(ind - 19) = y;

    v = (vdata(ind, 2) + vdata(ind + 1, 2)) * 0.25;
    dt = vdata(ind + 1, 1) - vdata(ind, 1);
    x = x + dt * v * cosd(81 - (vdata(ind, 4) + vdata(ind + 1, 4)) * 0.5);
    y = y + dt * v * sind(81 - (vdata(ind, 4) + vdata(ind + 1, 4)) * 0.5);
end

x2 = x2(:);
y2 = y2(:);
X2 = [ones(size(x2)) x2];
b2 = regress(y2, X2);

% WRITE ABOVE INTO FILE
for ind = 1:75
    if ind < 36
        writex = x1(ind);
        writey = x1(ind) * b1(2) + b1(1);
    else
        writex = x2(ind - 30);
        writey = x2(ind - 30) * b2(2) + b2(1);
    end
    
    fprintf(vlf, '%e   ', vdata(ind, 1)); % time
    fprintf(vlf, '%e   ', writex); % x pos
    fprintf(vlf, '%e   \n', writey); % y pos
end

% FINAL DIRECTION CHANGE
x0 = x;
y0 = y;
r = 0.01;
index = 76;
for time = 1:5:80
    theta = 83 - 52 / 80 * time;
    x = r * cosd(theta) + x0;
    y = r * sind(theta) + y0;
    
    fprintf(vlf, '%e   ', vdata(index, 1));
    fprintf(vlf, '%e   ', x);
    fprintf(vlf, '%e   \n', y);
    
    index = index + 1;
end

% EXTRA RUNTIME
fprintf(vlf, '%e   ', vdata(104, 1));
fprintf(vlf, '%e   ', x);
fprintf(vlf, '%e   \n', y);

fclose(vlf);


