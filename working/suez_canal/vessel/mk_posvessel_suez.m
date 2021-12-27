% PRODUCE VESSEL POSITION FILE FOR suez_canal

% - CREDIT -------------------------------------------------------------- 
% Maritime Casuality Specialist. (2021). Ever Given grounding in 
% Suez Canal - AIS based Dynamic Reconstruction. Retrieved July
% 7, 2021 at vimeo.com/531626438
%
% VesselFinder. (2021). EVER GIVEN, Container Ship - Details and current 
% position - IMO 9811000 MMSI 353136000 - VesselFinder. Retrieved July
% 7, 2021 at www.vesselfinder.com/vessels/EVER-GIVEN-IMO-9811000-MMSI-353136000

% - DESCRIPTION----------------------------------------------------------
% The following file will produce positional data of the Ever Given.
% The positional data will be generated from vessel_data.csv. 

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

for ind = 1:74
    v = (vdata(ind, 2) + vdata(ind + 1, 2)) * 0.25;
    dt = vdata(ind + 1, 1) - vdata(ind, 1);
    x = x + dt * v * cosd(81 - (vdata(ind, 4) + vdata(ind + 1, 4)) * 0.5);
    y = y + dt * v * sind(81 - (vdata(ind, 4) + vdata(ind + 1, 4)) * 0.5);

    fprintf(vlf, '%e   ', vdata(ind, 1)); % time
    fprintf(vlf, '%e   ', x); % x pos
    fprintf(vlf, '%e   \n', y); % y pos
end

% EXTRA RUNTIME
v = (vdata(75, 2) + vdata(76, 2)) * 0.25;
dt = vdata(75, 1) - vdata(76, 1);
x = x + dt * v * cosd(81 - (vdata(75, 4) + vdata(76, 4)) * 0.5);
y = y + dt * v * sind(81 - (vdata(75, 4) + vdata(76, 4)) * 0.5);

fprintf(vlf, '%e   ', vdata(104, 1)); % time
fprintf(vlf, '%e   ', x); % x pos
fprintf(vlf, '%e   \n', y); % y pos

% FINAL DIRECTION CHANGE
% x0 = x;
% y0 = y;
% r = 5;
% count = 0;
% for time = 1:20
%     theta = 40 + time / 20;
%   x = r * cosd(theta) + x0;
%    y = r * sind(theta) + y0;
%    
%    fprintf(vlf, '%e   ', 76 + count);
%    fprintf(vlf, '%e   ', x);
%    fprintf(vlf, '%e   \n', y);
%    
%    count = count + 1;
% end

fclose(vlf);


