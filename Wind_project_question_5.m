clc; clear;

% ---- Points 1–9 data ----
names = {'1: 10 kV','2: 10 kV bus','3: 20/10 kV side (20 kV)','4: 20 kV bus', ...
         '5: 70/20 kV side (70 kV)','6: 200/70 kV side (70 kV)', ...
         '7: 200/70 kV side (200 kV)','8: 200 kV bus','9: Strong grid'};
U   = [10 10 20 20 70 70 200 200 200]*1e3;             % V (L-L)
Sk  = [500 900 1100 1300 1300 1800 2200 3000 inf]*1e6; % VA
Ilo = [150 150 130 130 130 120  42  84  64];           % A
Ihi = [370 390 345 345 345 325 114 228 256];           % A
n   = [  2   4   2   2   2   1   2   1   1];           % number of lines
pfL = 0.9;  tanPhiL = tan(acos(pfL));                   % load pf

% PF range for the wind farm (unity included)
pfGrid = linspace(0.95,1.00,11);                        % 0.95 … 1.00
signs  = [-1,+1];                                       % -1=capacitive (leading), +1=inductive (lagging)

% Storage
P_Imin_pf1  = zeros(9,1);   P_Imax_pf1  = zeros(9,1);   P_cons_pf1  = zeros(9,1);
P_Imin_best = zeros(9,1);   P_Imax_best = zeros(9,1);   P_cons_best = zeros(9,1);
PF_Imin_used = zeros(9,1);  PF_Imax_used = zeros(9,1);
Mode_Imin = strings(9,1);   Mode_Imax = strings(9,1);
Limiting = strings(9,1);

for k = 1:9
    Un = U(k);
    Xth = 0; if ~isinf(Sk(k)), Xth = (Un^2)/Sk(k); end        % purely inductive per brief
    if Un >= 70e3, Vmin = 0.90*Un; Vmax = 1.10*Un;
    else,         Vmin = 0.90*Un; Vmax = 1.06*Un; end

    % Load powers at Imin/Imax
    Smin = n(k)*sqrt(3)*Un*Ilo(k);   PminL = pfL*Smin;   QminL = PminL*tanPhiL;
    Smax = n(k)*sqrt(3)*Un*Ihi(k);   PmaxL = pfL*Smax;   QmaxL = PmaxL*tanPhiL;

    % ---- Baseline: pf = 1.00 (Qw = 0) ----
    P_Imin_pf1(k) = Pmax_at_load(Un,Xth,PminL,QminL,Vmin,Vmax,0);
    P_Imax_pf1(k) = Pmax_at_load(Un,Xth,PmaxL,QmaxL,Vmin,Vmax,0);
    P_cons_pf1(k) = min(P_Imin_pf1(k), P_Imax_pf1(k));

    % ---- Best within pf ∈ [0.95, 1.00], capacitive & inductive ----
    [P_Imin_best(k), PF_Imin_used(k), Mode_Imin(k)] = best_at_load(Un,Xth,PminL,QminL,Vmin,Vmax,pfGrid,signs);
    [P_Imax_best(k), PF_Imax_used(k), Mode_Imax(k)] = best_at_load(Un,Xth,PmaxL,QmaxL,Vmin,Vmax,pfGrid,signs);
    P_cons_best(k) = min(P_Imin_best(k), P_Imax_best(k));
    Limiting(k) = string( ifelse(P_Imin_best(k) <= P_Imax_best(k), 'Imin', 'Imax') );
end

T_best = table(string(names)', U(:)/1e3, ...
    P_Imin_best(:), PF_Imin_used(:), Mode_Imin, ...
    P_Imax_best(:), PF_Imax_used(:), Mode_Imax, ...
    P_cons_best(:), Limiting, ...
    'VariableNames',{'Point','U_kV', ...
    'Pmax_Imin_MW','PF_Imin_used','Mode_Imin', ...
    'Pmax_Imax_MW','PF_Imax_used','Mode_Imax', ...
    'Pmax_cons_MW','Limiting_case'});

disp('--- Best within pf ∈ [0.95, 1.00], capacitive/inductive ---'); disp(T_best);

% ===== Helpers =====
function [Pbest,pfbest,modebest] = best_at_load(Un,Xth,Pload,Qload,Vmin,Vmax,pfGrid,signs)
    Pbest = 0; pfbest = 1.00; modebest = "unity";
    for pfw = pfGrid
        if abs(pfw-1.0) < 1e-12
            P = Pmax_at_load(Un,Xth,Pload,Qload,Vmin,Vmax,0);  % unity (Qw=0)
            if P > Pbest, Pbest=P; pfbest=1.00; modebest="unity"; end
        else
            tanPhiW = tan(acos(pfw));
            for s = signs
                P = Pmax_at_load(Un,Xth,Pload,Qload,Vmin,Vmax, s*tanPhiW);
                if P > Pbest
                    Pbest = P; pfbest = pfw;
                    if s < 0
                        modebest = "capacitive";
                    else
                        modebest = "inductive";
                    end
                end
            end
        end
    end
end

function Pmw = Pmax_at_load(Un,Xth,Pload,Qload,Vmin,Vmax,tanPhiW_sign)
    lo = 0; hi = 10;
    while true
        if ~is_feasible(Un,Xth,Pload,Qload,hi,tanPhiW_sign,Vmin,Vmax) || hi>2e4, break; end
        hi = hi*2;
    end
    while (hi - lo) > 0.1
        mid = 0.5*(hi+lo);
        if is_feasible(Un,Xth,Pload,Qload,mid,tanPhiW_sign,Vmin,Vmax), lo = mid; else, hi = mid; end
    end
    Pmw = lo;
end

function ok = is_feasible(Un,Xth,Pload,Qload,Pmw,tanPhiW_sign,Vmin,Vmax)
    Pw = Pmw*1e6;                 % W
    Qw = tanPhiW_sign*Pw;         % + = lagging (inductive), - = leading (capacitive)
    Snet = (Pload - Pw) + 1j*(Qload - Qw);
    V = Un + 0j;
    for it = 1:40
        I = conj(Snet)/(sqrt(3)*V);
        V = Un - 1j*Xth*I;
    end
    Vm = abs(V);
    ok = (Vm>=Vmin && Vm<=Vmax);
end

% Tiny inline helper for Limiting case label
function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end
