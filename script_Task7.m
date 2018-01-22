%%%% TASK 7 SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measuring the variation in the accuracy of the simulation performed for
% Task 6 as a function of the separation distance between the point
% sources. The accuracy is computed using the normalised rms error E (%).
%
% To compute a reference solution numerically by distributing a large
% number of point sources along the line source, at a source separation
% distance of 1 micron.
%
% Parameters from Task 6 are used.
%
% Pressure field evaluated on a grid defined as:
% Spatial grid:             -2mm <= x <= 2mm
%                         50 micron <= y <= 2mm
%                                   z = 0
% Spatial step size:     50 micron (in every direction)
% Temporal range:           0 <= t <= 3 microsecond
% Temporal step size:               10 ns
%
%
% Output movie:
% Displaying the pressure field at t=1 microsecond, each frame deplaying
% the pressure field obtained using a different point source separation
% distance. (ranging from 5 micron to 200 micron, in steps of 5 micron.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clearvars; % to avoid errors due to different parameters set in this task

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set initial parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set spatial grid [mm]
x = (-2:50E-3:2);
y = (50E-3:50E-3:2);
z = 0;
% set temporal axis and step size [s]
delta_t = 10E-9;
t = (0:delta_t:3E-6);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute reference solution %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% distribute large number of point sources along line source, with source
% separation distance of 1 micron
% set coordinates of point sources [mm]
xs_r = (-0.5:1E-3:0.5);
ys_r= 0;
zs_r = 0;

% return 1D array of all point source coordinates in the order of
% (xs,ys,zs)
ptsource_coords_r = reshape(combvec(xs_r,ys_r,zs_r),1,[])';

% preallocate size of total pressure field
tot_p_r = zeros(numel(x),numel(y));

% loop over all point source coordinates and compute each corresponding
% individual pressure field, then sum up all fields into tot_p
for n = 1:3:numel(ptsource_coords_r)
    
   % compute pressure field for each individual point source
   each_p_r = comp_press_field_point_source(c,p_0,x,y,z,ptsource_coords_r(n),ptsource_coords_r(n+1),ptsource_coords_r(n+2),t,delta_t);

   % sum all individual pressure fields
   tot_p_r = tot_p_r + each_p_r;
   
   % normalise total pressure field by dividing it by the number of point
   % sources used
   norm_tot_p_r = tot_p_r ./ (numel(ptsource_coords_r)/3);
   
end

% calculate excitation function
exc_fn = comp_Gaussian_tone_burst(f_0,sigma,delta_t);

% time index when t = 1 microsecond
ind_tspecific = find(t == 1E-6);

% preallocate size of total s (601x601)
tot_s_r = zeros(numel(x),numel(y));

% loop over all x and y values, compute each s corresponding to each grid
% point over all time points, sum all individual s at t=1 microsecond to
% tot_s
for ind_x = 1:numel(x)
    for ind_y = 1:numel(y)
        
        % get each pressure field at each grid point over all time points
        each_press_r = norm_tot_p_r(ind_x,ind_y,1,1:numel(t));
        
        % compute convolution between each pressure field and excitation
        % function over time
        each_s_r = conv(each_press_r(:),exc_fn,'same');
        
        % store each s at t=1 microsecond into tot_s, which gives a 601x601
        % array
        tot_s_r(ind_x,ind_y) = each_s_r(ind_tspecific,:);
        
        
    end
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine point source separation distances & coordinates %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define all sets of point source separation distances
% range of distances step size [mm]
delta_dist = 5E-3;
% range of distances [mm]
ptsource_dist = (5E-3:delta_dist:200E-3);


% determine xs coordinates [mm]%

% find the maximum number of xs coordinates corresponding to the smallest
% point source separation distance
max_xs_no = numel(-0.5:5E-3:0.5);

% preallocate cell size for all xs coordinates 
xs = cell(numel(ptsource_dist),max_xs_no);


% set ys and zs coordinates [mm]
ys = 0;
zs = 0;


% loop over all possible point source separation distances, obtain a cell
% containing sets of xs coordinates for all different separation distances
for n = 1:numel(ptsource_dist)
   
    % store each set of xs coordinates into a cell
    xs{n} = (-0.5:ptsource_dist(n):0.5);
 
end


% find the maximum number of all possible xs,ys,zs values, which corresponds
% to the smallest point source separation distance
max_coords_no = numel(reshape(combvec(xs{1},ys,zs),1,[])');

% preallocate cell size for all point source coordinates
ptsource_coords = cell(numel(ptsource_dist),max_coords_no);

% loop over all point source separation distances, obtain a cell containing
% sets of point source coordinates for all different separation distances
for n = 1:numel(ptsource_dist)
    % return 1D array of all point source coordinates in the order of
    % (xs,ys,zs)
    ptsource_coords{n} = reshape(combvec(xs{n},ys,zs),1,[])';
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute pressure fields for different separation distances %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open one figure to display all frames one at a time
figure;


% open videowriter, name video and set format
vid = VideoWriter('task7.mp4','MPEG-4');
vid.FrameRate = 8; % set frame rate to 8 frames per second
open(vid);



%%%% Plot reference solution first %%%%
% note: plotting the transpose of tot_s
% 2D filled contour plot with 200 contour levels set
contourf(x,y,tot_s_r',200,'edgecolor','none')
        
title(['Source spacing: ' num2str(1) '\mum' ', reference solution']) % update title for each frame
xlabel('x [mm]')
ylabel('y [mm]')
axis ij % reverse y axis
axis equal % use equal unit lengths along both axes
colormap(hot)
colorbar
caxis([-30 30]) % set colorbar axis limits (corresponding to pressure values)
        
% get current frame and write this frame into the video
frame = getframe(gcf);
writeVideo(vid,frame);



%%%% For all other different point source separation distances %%%%

% preallocate sizes of tot_p and tot_s
tot_p = zeros(numel(x),numel(y));
tot_s = zeros(numel(x),numel(y));


% loop over each point source separation distance, display each pressure
% field corresponded in each frame
for n = 1:numel(ptsource_dist)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute p(x,y,z,t) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % loop over all possible point source coordinates to compute pressure
    % field
    for m = 1:3:numel(ptsource_coords{n})
    
        % compute pressure field for each individual point source
        each_p = comp_press_field_point_source(c,p_0,x,y,z,ptsource_coords{n}(m),ptsource_coords{n}(m+1),ptsource_coords{n}(m+2),t,delta_t);

        % sum all individual pressure fields
        tot_p = tot_p + each_p;
        
        % normalised total pressure field by dividing it by the number of
        % point sources used
        norm_tot_p = tot_p./(numel(ptsource_coords{n})/3);
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compute s(x,y,z,t) %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % loop over all x and y values to compute s
    for ind_x = 1:numel(x)
        for ind_y = 1:numel(y)
        
            % get each pressure field at each grid point over all time points
            each_press = norm_tot_p(ind_x,ind_y,1,1:numel(t));
        
            % compute convolution between each pressure field and excitation
            % function over time
            each_s = conv(each_press(:),exc_fn,'same');
        
            % store each s at t=1 microsecond into tot_s, which gives a 601x601
            % array
            tot_s(ind_x,ind_y) = each_s(ind_tspecific,:);
        
        end
    end
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Calculate normalied RMS error %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % expectation value of reference solution s squared
    exp_s_ref_squared = mean2((tot_s_r).^2);
    
    % expectation value of each current norm_tot_s (corresponding to each source
    % separation distance) minus exp_s_ref, all squared
    exp_s_diff_squared = mean2((tot_s - tot_s_r).^2);
    
    % calculate normalised RMS error E
    norm_RMS_err = 100 * (sqrt(exp_s_diff_squared) / sqrt(exp_s_ref_squared));
    
        
        
    %%%%%%%%%%%%%%%%%%
    %%%% Plotting %%%%
    %%%%%%%%%%%%%%%%%%
        
    % note: plotting the transpose of tot_s
    % 2D filled contour plot with 200 contour levels set
    contourf(x,y,tot_s',200,'edgecolor','none')
        
    % update title for each frame, displaying source spacing in um and
    % error to 2 decimal places
    title(['Source spacing: ' num2str(ptsource_dist(n)*10^3) '\mum' ', error = ' sprintf('%.2f',norm_RMS_err) '%'])
    xlabel('x [mm]')
    ylabel('y [mm]')
    axis ij % reverse y axis
    axis equal % use equal unit lengths along both axes
    colormap(hot)
    colorbar
    caxis([-30 30]) % set colorbar axis limits (corresponding to pressure values)
        
    % get current frame and write this frame into the video
    frame = getframe(gcf);
    writeVideo(vid,frame);
        
    % reset total pressure field to zero for the next separation distance
    tot_p = zeros(numel(x),numel(y));
    tot_s = zeros(numel(x),numel(y));
    
end

close(vid); % close videowriter after loop ends









