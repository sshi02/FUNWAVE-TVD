% MAKE MOVIE FOR suez_canal

% - DESCRIPTION---------------------------------------------------------- 
% The following script will produce a .avi movie for suez_canal, using the
% given data provided by mk_depth_suez and mk_vessel_suez. 

clear all; 
vdir = '../vessel/';
fdir = '../output/';

% figure setup
vdata = csvread([vdir 'vessel_data.csv'], 1, 0);
eta = load([fdir 'eta_00001']);

[n, m] = size(eta);
dx = 3.0;
dy = 3.0;
x = [0:m - 1] * dx;
y = [0:n - 1] * dy;

nfile = 0:1:599; % 600 data files

colormap jet;
wid = 1;
len = 18;
set(gcf, 'units', 'inches', 'paperunits', 'inches', 'papersize', ...
    [wid len], 'position', [1 1 wid len], 'paperposition', [0 0 wid len]);
clf;

vidObj = VideoWriter('suez_mov.avi');
vidObj.FrameRate = 12;
open(vidObj);

for num = 1:2:length(nfile)
    fnum = sprintf('%.5d', nfile(num));
    eta = load([fdir 'eta_' fnum]);
   
    % NOTE ---- REMOVE WHEN DONE ----------------------------------------
    % Likely should map on satillite data first, then plot the surface or
    % psuedo-color map. Unclear whether or not using surf or pcolor is
    % better orientation, clarity, and looks-wise. Will have to tweak caxis
    % and orientation a lot. 
    
    pcolor(x, y, eta), shading flat
    hold on;
    caxis([-2.1 2.1]);
    title([' Time = ' num2str(num) ' sec ']);
    xlabel(' x (m) ');
    ylabel(' y (m) ');
    
    f = getframe(gcf);
    writeVideo(vidObj, f);
end

close(vidObj);
    