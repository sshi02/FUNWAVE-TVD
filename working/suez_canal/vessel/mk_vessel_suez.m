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
% The positional data will be generated from mk_posvessel_suez. Due to 
% complications between dense data and FUNWAVE, the positional data will 
% be smoothed out.

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

% first major path
pos1 = [0.000000e+00   1.985318e+02   2.529673e+02; ... $ intial pos
        2.360000e+02   9.723960e+01   1.788378e+03; ... $ start dir change
        2.410000e+02   9.681746e+01   1.820626e+03; ...   
        2.460000e+02   9.653603e+01   1.852874e+03; ...
        2.510000e+02   9.647974e+01   1.885124e+03; ... $ xvelocity flip
        2.560000e+02   9.653603e+01   1.917374e+03; ...   
        2.610000e+02   9.667729e+01   1.949749e+03; ...   
        2.660000e+02   9.696090e+01   1.982248e+03; ... $ end dir change   
        ];
    
% second major path
pos2 = [2.870000e+02   1.000155e+02   2.118305e+03; ... $ small slow
        3.790000e+02   1.899764e+02   2.613157e+03];
        

    % WRITE ABOVE
for ind = 1:size(pos1, 1)
    fprintf(vlf, '%e   %e   %e\n', pos1(ind, :));
end

for ind = 1:size(pos2, 1)
    fprintf(vlf, '%e   %e   %e\n', pos2(ind, :)); 
end

% FINAL DIRECTION CHANGE
pos3 = zeros(16);
x0 = 1.899764e+02;
y0 = 2.613157e+03;
r = 5;
index = 76;
for time = 1:20
    theta = 173 - 7 / 80 * time;
    x = r * cosd(theta) + x0 + 5;
    y = r * sind(theta) + y0;
    
    fprintf(vlf, '%e   ', vdata(index, 1));
    fprintf(vlf, '%e   ', x);
    fprintf(vlf, '%e   \n', y);
    
    index = index + 1;
end


% EXTRA RUNTIME
fprintf(vlf, '%e   ', vdata(104, 1));
fprintf(vlf, '%e   ', x + 1);
fprintf(vlf, '%e   \n', y + cos(3.14 / 6));

fclose(vlf);


