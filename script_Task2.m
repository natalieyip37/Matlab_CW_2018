%%%% TASK 2 SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using function comp_Gaussian_tone_burst to compute a Gaussian tone burst
% for the following parameters:
%
% f_0 = 10 MHz
% sigma = 50 ns
% delta_t  10 ns
%
% Using approximation to the Gaussian: -4*stdev <= t <= +4*stdev
%
%
% Return plot - amplitude of Gaussian tone burst against time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clearvars;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set given initial parameters %%%%
% NOTE: SI units for all variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set centre frequency [Hz]
f_0 = 10E6;
% set standard deviation [s]
sigma = 50E-9;
% set temporal step size [s]
delta_t = 10E-9;
% set time range using approximation -4*stdev <= t <= +4*stdev
t = (-4*sigma:delta_t:4*sigma);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Use function to calculate Gaussian tone burst %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gauss_tb = comp_Gaussian_tone_burst(f_0,sigma,delta_t);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot amplitude of Gaussian tone burst against time %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(t,gauss_tb)
xlabel('Time [s]')
ylabel('Amplitude [a.u.]')
title('Normalised Gaussian tone burst')


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function checking %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate integral of Gaussian tone burst over time range
syms t
int_gauss_tb = int(exp(-t.^2 ./ (2*sigma.^2)) .* sin(2*pi.*f_0.*t),t(1),t(end));

% check if integral of Gaussian tone burst gives zero (since it is
% antisymmetric)
if int_gauss_tb == 0
    disp('Output is correct since integral of Gaussian tone burst gives zero')
else
    disp('Output is incorrect since integral of Gaussian tone burst does not give zero')
end








