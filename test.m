load('test_data.mat')

%%% Enter you own data here: %%%
% D := rotor diameter in meter
% k := wake decay constant, set to k=0.0750
% u := wind speeds in m/s
% dir := wind directions in degree
% Ct := turbine thrust coefficient
% Ct_u := wind speed vector
% xy := geographical alignment of wind turbines
%       xy(:,1) := x-axis coordinates
%       xy(:,2) := y-axis coordinates
%       size(xy,1) = number of turbines

D =  ;              % rotor diameter in meter
k = 0.0750;         % wake decay constant, set to k=0.0750
u =  ;              % wind speeds in m/s
dir =    ;          % wind directions in degree
Ct =   ;            % turbine thrust coefficient
Ct_u =     ;        % wind speed vector
xy = [   ,   ];     % geographical alignment of the ? wind turbines
      

%% Run this code to get the reduced wind speeds at each turbine

u2 = reduced_wind_speeds(xy,D,k,u,dir,Ct,Ct_u);
subplot(1,2,1)
scatter(xy(:,1),xy(:,2))
xlabel('W-E coordinate (m)')
ylabel('S-N coordinate (m)')
subplot(1,2,2)
scatter(dir,u-mean(u2,2))
xlabel('Direction (degrees)')
ylabel('Reduction in wind speed (m/s)')