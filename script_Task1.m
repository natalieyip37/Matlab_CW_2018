%%%% TASK 1 SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Using function comp_press_field_point_source to compute the pressure at
% grid point (x,y,z) = (1,1,0) mm for all time sample 0<=t<=3
% microsecond using a temporal step size of delta_t = 10 nanosecond.
%
% Initial parameters given:
% xs = ys = zs = 0 mm
% c = 1500 m/s (simulating water or soft biological tissue)
% p_0 = 1 Pa*m
%
% Return plot - resulting pressure [Pa] against time [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

clearvars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set given initial parameters %%%%
%%%  NOTE units for all variables  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set pressure grid point [mm]
x = 1; y = 1; z = 0;
% set point source grid point [mm] 
xs = 0; ys = 0; zs = 0;
% set speed of sound [mm/s]
c = 1500E3;
% set initial acoustic pressure [Pa*mm]
p_0 = 1E3;
% set temporal step size [s]
delta_t = 10E-9;
% create array of time points in range 0<=t<=3 with temporal step size
% delta_t [s]
t = (0:delta_t:3E-6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Use function to calculate pressure field (giving 4D array) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1_4D = comp_press_field_point_source(c,p_0,x,y,z,xs,ys,zs,t,delta_t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot resulting pressure against time %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(t,p1_4D(:))
xlabel('Time [s]')
ylabel('Pressure [Pa]')
title('Pressure as function of time for an acoustic point source')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function checking %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate pressure value using the pressure equation given
r = sqrt((x-xs)^2 + (y-ys)^2 + (z-zs)^2);
cal_p = p_0 / (4*pi*r);

% check if function gives the same pressure value as equation
if abs(cal_p - max(p1_4D(:))) <= eps(cal_p)
    disp('Function and equation give the same pressure value, function is working correctly')
else
    disp('Different pressure values obtained from function and equation, please check function')
end


% calculate arrival time of acoustic pulse using equation given
cal_t = r / c;
round_cal_t = round(cal_t,2,'significant'); % round calculated arrival time to 2 sig fig

% find arrival time from function
arr_t_ind = find(p1_4D == max(p1_4D(:)));
arr_t = t(arr_t_ind);

% check if function gives the same arrival time as the equation
if abs(round_cal_t - arr_t) <= eps(round_cal_t)
    disp('Function and equation give the same arrival time of acoustic pulse, function is working correctly')
else
    disp('Different arrival times obtained from function and equation, please check function')
end






