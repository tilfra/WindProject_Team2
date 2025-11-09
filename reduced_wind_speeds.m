function [u2] = reduced_wind_speeds(xy,D,k,u,dir,Ct,Ct_u)
% REDUCED_WIND_SPEEDS - Jensen-like wake model (interp1 for Ct, no fit)
% Inputs:
%   xy    : num x 2 coordinates of turbines
%   D     : rotor diameter (m)
%   k     : wake decay constant
%   u     : time x num matrix (or Nx1 -> replicated)
%   dir   : time x 1 wind direction degrees (0..359)
%   Ct    : Ct vector (same length as Ct_u)
%   Ct_u  : speeds vector for Ct
%
% Output:
%   u2    : reduced wind speeds (time x num)


% Ct interpolation (no toolbox). clamp to [0,1] for safety.
Ct_fun = @(v) max(0, min(1, interp1(Ct_u, Ct, v, 'linear', 'extrap')));

num = size(xy,1); % Number of WTs
d = pdist2(xy,xy); % distance between WTs
dir = round(dir); 
dir(dir==0) = 360;
if size(u,2) == 1
    u = repmat(u,1,num); % same raw wind speed for all WTs
end

% Which WTs are within wake for different angles?
pin = false(num,360,num); % initialize logical 
for n = 1:num
    pin = in_wake(xy,n,k,D,pin); % within wake for angles 1-360
end

% Calculate catic = 1 - V/U (use Ct_fun)
catic = zeros(num,size(u,1),num,'single');
for n = 1:num   
    for m = 1:num
        if n~=m
            angles = find(squeeze(pin(m,:,n)));    % upwind angles where m affects n
            if isempty(angles), continue; end
            temp = ismember(dir, angles);         % times when wind comes from those angles
            if ~any(temp), continue; end
            % use Ct_fun (interpolation) for the speeds at those times
            cu = Ct_fun(u(temp,n));               % vector
            % Jensen-style expression (original form preserved)
            catic(m,temp,n) = (1 - sqrt(1 - cu)) ./ (1 + 2*k*d(n,m)/D).^2;
        end
    end
end

% Calculate reduced wind speeds
u2 = zeros(size(u),'like',u);
for n = 1:num
    vu = 1 - sqrt(sum(catic(n,:,:).^2,3));  % 1 x time
    vu = real(vu);
    vu(vu < 0) = 0;
    u2(:,n) = u(:,n) .* vu';
end

end
