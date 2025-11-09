% Wind directions in degrees
theta = 0:30:330;  

% Relative frequencies
f = [0.06 0.04 0.04 0.05 0.05 0.06 0.14 0.21 0.09 0.05 0.07 0.12];  

% Convert degrees to radians for trigonometric functions
theta_rad = deg2rad(theta);

% Compute weighted sums of sine and cosine components
x = sum(f .* cos(theta_rad));
y = sum(f .* sin(theta_rad));

% Compute circular mean (in radians)
theta_mean_rad = atan2(y, x);

% Convert back to degrees
theta_mean_deg = rad2deg(theta_mean_rad);

% Ensure result is between 0 and 360°
if theta_mean_deg < 0
    theta_mean_deg = theta_mean_deg + 360;
end

disp(['Mean wind direction: ', num2str(theta_mean_deg), '°']);