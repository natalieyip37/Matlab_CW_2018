%%%% TASK 6 SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extend similation to acoustic sources of finite dimensions, using the
% Huygen's principle.
% 
% To model an acoustic line source located at (|x|<=0.5mm, y=z=0) as a
% finite collection of point sources.
% Each of these point sources coincide with a grid point.
%
% Parameters set:
% sound speed:                 1500 m/s
% initial pressure:            1 Pa*m
% spatial grid:            -2mm <= x <= 2mm
%                           0mm <= y <= 4mm
%                                 z = 0
% spatial step size:     50 micron (in every direction)
% temporal range:          0 <= t <= 4 microseconds
% temporal step size:             10 ns
%
% Gaussian tone burst parameters as set in Task 2:
% f_0 = 10 MHz
% sigma = 50 ns
% delta_t = 10 ns
%
%
% OUTPUT FIGURE:
% Display total pressure field observed at t = 1 microsecond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clearvars; % to avoid errors due to different parameters set for this task

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set initial parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set spatial grid [mm]
x = (-2:50E-3:2);
y = (0:50E-3:4);
z = 0;
% indices for x,y,z coordinates
ind_x = 1:numel(x);
ind_y = 1:numel(y);
ind_z = 1:numel(z);
% set temporal axis and step size [s]
delta_t = 10E-9;
t = (0:delta_t:4E-6);
ind_t = 1:numel(t);
% set speed of sound [mm/s]
c = 1500E3;
% set initial pressure [Pa*mm]
p_0 = 1E3;


% set Gaussian tone burst parameters from Task 2

% set centre frequency [Hz]
f_0 = 10E6;
% set standard deviation [s]
sigma = 50E-9;
% set temporal step size [s]
delta_t = 10E-9;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Generate pressure field %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine all grid points within line source (point source coordinates
% xs,ys,zs)
xs = x(abs(x) <= 0.5);
ys = y(y == 0);
zs = z(z == 0);

% return 1D array of all point source coordinates in the order of
% (xs,ys,zs)
ptsource_coords = reshape(combvec(xs,ys,zs),1,[])';

% preallocate size of total pressure field
tot_p = zeros(numel(x),numel(y));

% loop over all point source coordinates and compute each corresponding
% individual pressure field, then sum up all fields into tot_p
for n = 1:3:numel(ptsource_coords)
    
   % compute pressure field for each individual point source
   each_p = comp_press_field_point_source(c,p_0,x,y,z,ptsource_coords(n),ptsource_coords(n+1),ptsource_coords(n+2),t,delta_t);

   % sum all individual pressure fields
   tot_p = tot_p + each_p;
   
end

% calculate excitation function
exc_fn = comp_Gaussian_tone_burst(f_0,sigma,delta_t);

% time index when t = 1 microsecond
ind_tspecific = find(t == 1E-6);

% preallocate size of total s (601x601)
tot_s = zeros(numel(x),numel(y));

% loop over all x and y values, compute each s corresponding to each grid
% point over all time points, sum all individual s at t=1 microsecond to
% tot_s
for ind_x = 1:numel(x)
    for ind_y = 1:numel(y)
        
        % get each pressure field at each grid point over all time points
        each_press = tot_p(ind_x,ind_y,1,1:numel(t));
        
        % compute convolution between each pressure field and excitation
        % function over time
        each_s = conv(each_press(:),exc_fn,'same');
        
        % store each s at t=1 microsecond into tot_s, which gives a 601x601
        % array
        tot_s(ind_x,ind_y) = each_s(ind_tspecific,:);
        
        
    end
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Displaying total pressure field %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on

% note: plotting the transpose of tot_s
% 2D filled contour plot with 200 contour levels set
contourf(x,y,tot_s',200,'edgecolor','none')

title('Total pressure field [Pa]')
xlabel('x [mm]')
ylabel('y [mm]')
axis ij % reverse y axis
axis equal % use equal unit lengths along both axes
colormap(hot)
colorbar
caxis([-400 400]) % set colorbar limits (corresponding to pressure values)


