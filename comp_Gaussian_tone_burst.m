function norm_gauss_tb = comp_Gaussian_tone_burst(f_0,sigma,delta_t)
% gauss_tb = comp_Gaussian_tone_burst(f_0,sigma,delta_t)
%
% This function computes a normalised Gaussian tone burst for arbitrary input values
% of centre frequency, standard deviation, and temporal step size.
% Using approximation: -4*stdev <= t <= +4*stdev
%
% INPUTS:
%       f_0:            centre frequency [Hz]
%       sigma:          standard deviation of the Gaussian envelope
%                       (exp(-t^2/2*sigma^2)), gives the width in time of
%                       the tone burst [s]
%       delta_t:        temporal step size [s] 
%
% OUTPUT:
%       norm_gauss_tb:  1D array - amplitudes of normalised Gaussian tone
%                       burst [a.u.]
%
%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Error checking %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

for input = [f_0,sigma,delta_t]
    % check if all inputs are positive
    if input < 0
        error('All inputs have to be positive')
    end
    % check if all inputs are numeric, real and scalar
    if ~isnumeric(input) || ~isreal(input) || ~isscalar(input)
        error('All inputs have to be numeric, real and scalar')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create time range %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create an array of time points using approximation -4*stdev <= t <=
% +4*stdev, and temporal step size (delta_t)
t = (-4*sigma:delta_t:4*sigma);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate amplitude of normalised Gaussian tone burst over t %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% output 1D array of Gaussian tone burst amplitudes over all time points
norm_gauss_tb = exp(-t.^2 ./ (2*sigma.^2)) .* sin(2*pi.*f_0.*t);






