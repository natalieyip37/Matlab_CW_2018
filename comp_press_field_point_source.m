function p_4D = comp_press_field_point_source(c,p_0,x,y,z,xs,ys,zs,t,delta_t)
% p_4D = comp_press_field_point_source(c,p_0,x,y,z,xs,ys,zs,t,delta_t)
%
% 
% This function computes the pressure field generated by an acoustic point
% source as a function of time over a 3D grid of sample points.
%
% INPUTS:
%       c:              speed of sound [mm/s]
%       p_0:            initial pressure amplitude [Pa*mm] (*Note unit*)
%       x,y,z:          coordinates (x,y,z) of grid point where
%                       pressure is to be evaluted [mm] 
%       xs,ys,zs:       cooridinates (x_s,y_s,z_s) of acoustic
%                       point source [mm]
%       t:              1D array - time points for which the field is computed [s]
%
% OUTPUTS:
%       p_4D:           4D array - pressure field (x,y,z,t)
%                       generated by an acoustic point
%                       source as a function of time over a 3D grid of
%                       sample points [Pa] 
%
%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error checking %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% check if all inputs are numeric and real
for input = [c,p_0,x,y,z,xs,ys,zs,t,delta_t]
   if ~isnumeric(input) || ~isreal(input)
       error('All inputs have to be numeric and real')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Obtain initial parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define indices for all time points
ind_t = 1:numel(t);

% coordinates of all grid points where pressure will be evaluated
% return a 1D vector
p_coords = reshape(combvec(x,y,z),1,[])';

% define indices for each grid dimension (x,y,z)
ind_x = 1:numel(x);
ind_y = 1:numel(y);
ind_z = 1:numel(z);

% obtain a 1D vector of all indices in the order of x,y,z
ind = reshape(combvec(ind_x,ind_y,ind_z),1,[])';

% calculate constant variables in pressure equation
const = p_0 / (4*pi);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate pressure field %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate size of pressure 4D array
p_4D = zeros(numel(x),numel(y),numel(z),numel(t));


% loop over all possible pressure coordinates (x,y,z)
for n = 1:3:numel(p_coords)
    
    % define each grid point 
    grid_pt = [p_coords(n),p_coords(n+1),p_coords(n+2)];
    
    % calculate distance between each grid point and point source
    r = sqrt((grid_pt(1) - xs)^2 + (grid_pt(2) - ys)^2 + (grid_pt(3) - zs)^2)';
    
    % dirac delta function, return a logical: 1 if t = r/c, otherwise 0
    dirac = bsxfun(@le,abs(t-r./c),delta_t/2);
    
    % if r is zero, pressure obtained from equation will give NaN,
    % this sets pressure to zero if r is zero (which corresponds to the
    % point of origin),
    % otherwise calculate pressure field using equation (for r~=0)
    if r == 0
        p = 0;
    else
        % calculate pressure field
        p = const * (dirac./r);
    end
    
    % store each pressure at each grid point, corresponding to each
    % dimension (x,y,z,t), returning a 4D double array.
    % for each grid point (x,y,z), all time points are evaluated
    p_4D(ind(n),ind(n+1),ind(n+2),ind_t) = p;
    
    

end




