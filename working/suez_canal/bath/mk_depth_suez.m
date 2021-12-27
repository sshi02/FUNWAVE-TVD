% PRODUCE DEPTH DATA FILE FOR suez_canal

% - CREDIT -------------------------------------------------------------- 
% Peterson, Curt & Butler, Virginia & Feathers, James & Cruikshank,
% Kenneth. (2014). Geologic Records of Net Littoral Drift, Beach Plain
% Development, and Paleotsunami Runup, North Sand Point, Olympic
% Washington, USA. Northwest Science. 88. 314-328. 10.3955/046.088.0406.
%
% Suez Canal Authority. Suez Canal Cross Sections. Retrieved July 7, 2021
% at web.archive.org/web/20150824031507/http://www.suezcanal.gov.eg/Files/Circular/Suez%20Canal%20Cross%20Section.pdf

% - DESCRIPTION---------------------------------------------------------- 
% The following script will produce a bathymetry file for suez_canal, using
% a diagram which shows the shape and size of the canal as follows:
% - Width Extrema: 121 m, 313 m
% - Depth Extrema: 0 m , 24 m
% - Slope: 4:1
% The desired dimensions in this file will be 104 x 1500, defined below.
% The script will write to file 'suez_depth'

clear all;

l = 1500; % length
w = 104; %width
wmin = 40;
sl = 0.78;
depmax = 25;

% produce a row
a = 0:w;
for ind = 1:w + 1
    if a(ind) == 0 || a(ind) == w
        a(ind) = 0;
    elseif a(ind) < (w - wmin) / 2
        a(ind) = a(ind) * sl;
    elseif a(ind) > w - ((w - wmin) / 2)
        a(ind) = (w - a(ind)) * sl;
    else
        a(ind) = depmax;
    end
end

% write to file
dpf = fopen('suez_depth', 'w');
c = 1;
for ind = 1:l * (w + 1)
    if c == (w + 1)
        fprintf(dpf, ' %e\n', a(c));
        c = 1;
    else
        fprintf(dpf, ' %e', a(c));
        c = c + 1; 
    end
end
fclose(dpf);
    