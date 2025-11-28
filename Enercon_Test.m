% Enercon E-138 (approx power curve via cubic model)
clear; clc; close all;

%% ---------------- User / farm parameters ----------------
D = 138;                % rotor diameter [m] (Enercon E-138)
k = 0.075;              % wake decay constant (same as before)
rated_kW = 4260;        % turbine rated power (Enercon ~4.26 MW)
farmMW = 80;            % target farm size (MW)
cutin = 2.5;            % cut-in (m/s)
cutout = 28.0;          % cut-out (m/s)
v_r = 10.8;             % rated wind speed (m/s)

%% ---------------- power curve (cubic approximation) ----------------
P_u = (0:0.5:30)';     
P_kW = zeros(size(P_u));

% cubic region: v_ci <= v < v_r
mask = (P_u >= cutin) & (P_u < v_r);
P_kW(mask) = rated_kW .* ((P_u(mask).^3 - cutin^3) ./ (v_r^3 - cutin^3));

% plateau region: v_r <= v < cutout
mask2 = (P_u >= v_r) & (P_u < cutout);
P_kW(mask2) = rated_kW;

% safety clipping
P_kW = max(0, min(P_kW, rated_kW));

% interpolation function (returns kW)
Pfun = @(x) max(0, interp1(P_u, P_kW, x, 'linear', 0));

%% ---------------- Load wind speeds ----------------
T = readtable('Windspeeds(1).csv','PreserveVariableNames',true,'Delimiter',';');
timeRaw = T{:,1};
try
    timeVec = datetime(string(timeRaw),'InputFormat','yyyy-MM-dd HH:mm');
catch
    timeVec = datetime(string(timeRaw)); % robust parse
end
raw = T{:,2};
if iscell(raw) || isstring(raw)
    raw = strrep(string(raw),',','.');
    u_all = double(str2double(raw));
else
    u_all = double(raw);
end

u_col = u_all;   % use all years


H = numel(u_col);
fprintf('Loaded %d hourly wind speeds (all years)\n', H);

%% ---------------- Wind rose (12-sector): it is used for sector weighted AEP ----------------
theta = (0:30:330)'; 
freq  = [0.08 0.04 0.04 0.05 0.05 0.06 0.14 0.21 0.09 0.05 0.07 0.12]';
freq = freq / sum(freq);

%% ---------------- load Ct data from test_data.mat ----------------
S = load('test_data.mat');   
Ct_u = double(S.Ct_u(:));
Ct   = double(S.Ct(:));

%% ---------------- build layouts ----------------
nT = round((farmMW*1000)/rated_kW);   % number turbines (Enercon)
spacing = 7 * D;                      % 7*D spacing, same approach as Vestas script

% Layout A: offset / staggered (two-row)
xpos = (0:nT-1)' * spacing;
ypos = mod(0:nT-1,2)' * (4*D);   
xy_offset = [xpos, ypos];

% Layout B: arced (aligned with mean wind direction)
mean_dir = 222.7085;               % from wind rose (same as Vestas script)
theta_deg_arc = 90;                %  60/90/120
row_orient = mean_dir - 90;   % ή +90, δοκίμασε και δες στο plot
xy_arced = create_arced_layout(nT, spacing, theta_deg_arc, row_orient, [0 0]);

layouts = {"offset", xy_offset; "arced", xy_arced};

%% ---------------- evaluate each layout (sector-weighted AEP) ----------------
results = table('Size',[2,4],'VariableTypes',{'string','double','double','double'}, ...
    'VariableNames',{'Layout','AEP_MWh_perYear','CapacityFactor_pct','WakeLoss_pct'});

for L=1:size(layouts,1)
    name = layouts{L,1};
    xy = layouts{L,2};
    nT_loc = size(xy,1);
    AEP_kWh = 0;
    for s = 1:numel(theta)
        dir_i = theta(s) * ones(H,1);   % constant direction for this sector
        u2 = reduced_wind_speeds(xy, D, k, u_col, dir_i, Ct, Ct_u);  % H x nT
        u2(u2 < cutin | u2 > cutout) = 0;
        Pk = Pfun(u2);
        AEP_i_kWh = sum(Pk,'all');  % kW-hours for the H timestamps in this sector
        AEP_kWh = AEP_kWh + freq(s) * AEP_i_kWh;
    end
    % Convert weighted sum to yearly MWh accounting for H (hours present)
    AEP_MWh = (AEP_kWh / 1000) * (8760 / H);
    P_single = Pfun(u_col);          
    AEP_no_wake_kWh = sum(P_single) * nT_loc;
    AEP_no_wake_MWh = (AEP_no_wake_kWh / 1000) * (8760 / H);
    rated_MW = rated_kW/1000;
    CF_pct = 100 * AEP_MWh / (rated_MW * nT_loc * 8760);
    CF_no_pct = 100 * AEP_no_wake_MWh / (rated_MW * nT_loc * 8760);
    WakeLoss_pct = 100 * (1 - AEP_MWh / AEP_no_wake_MWh);
    results.Layout(L) = string(name);
    results.AEP_MWh_perYear(L) = AEP_MWh;
    results.CapacityFactor_pct(L) = CF_pct;
    results.WakeLoss_pct(L) = WakeLoss_pct;
    fprintf('Layout %s: AEP = %.1f MWh/yr | CF = %.2f %% | wake loss = %.2f %%\n', ...
        name, AEP_MWh, CF_pct, WakeLoss_pct);
end

fprintf('Saved offset_vs_arced_comparison_enercon.csv\n');

%% ---------------- plot ----------------
figure('Name','Layouts (Enercon)'); 
subplot(1,2,1); scatter(xy_offset(:,1),xy_offset(:,2),70,'k','filled'); axis equal; title('Offset (staggered)'); xlabel('x (m)'); ylabel('y (m)');
subplot(1,2,2); scatter(xy_arced(:,1),xy_arced(:,2),70,'k','filled'); axis equal; title(sprintf('Arced (orient=%.2f°)', mean_dir)); xlabel('x (m)'); ylabel('y (m)');

%% ---------------- subfunction: create_arced_layout ----------------
function xy = create_arced_layout(n, spacing, theta_deg, orientation_deg, center)
% create N turbine coordinates along an arc
if nargin < 5 || isempty(center), center = [0,0]; end
if n < 2, error('n must be >= 2'); end
theta = deg2rad(theta_deg);
dphi = theta / (n-1);
r = spacing / dphi;
phi = linspace(-theta/2, theta/2, n);
x_local = r * sin(phi);
y_local = r * (1 - cos(phi));
rot = [cosd(orientation_deg) -sind(orientation_deg); sind(orientation_deg) cosd(orientation_deg)];
xy_rot = (rot * [x_local; y_local])';
xy = xy_rot + center(:)';
end
