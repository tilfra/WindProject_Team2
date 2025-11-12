clc
clear

names = {'1: 10 kV','2: 10 kV bus','3: 20/10 kV side (20 kV)','4: 20 kV bus', ...
         '5: 70/20 kV side (70 kV)','6: 200/70 kV side (70 kV)', ...
         '7: 200/70 kV side (200 kV)','8: 200 kV bus','9: Strong grid'};

U   = [10 10 20 20 70 70 200 200 200]*1e3;            % V (line-line)
Sk  = [500 900 1100 1300 1300 1800 2200 3000 inf]*1e6;% VA short-circuit
Ilo = [150 150 130 130 130 120  42  84  64];          % A (per line set)
Ihi = [370 390 345 345 345 325 114 228 256];          % A
n   = [  2   4   2   2   2   1   2   1   1];          % number of lines
pfL = 0.9; phi = acos(pfL);                           % load pf

Pmin = zeros(numel(U),1);   % MW at Imin
Pmax = zeros(numel(U),1);   % MW at Imax
Pcons = zeros(numel(U),1);  % MW conservative = min(Pmin,Pmax)

for k = 1:numel(U)
    Un = U(k);
    Xth = isinf(Sk(k)) * 0 + (~isinf(Sk(k))) * (Un^2 / Sk(k)); % ohms, inductive
    if U(k)/1e3 >= 70, Vmin = 0.90*Un; Vmax = 1.10*Un; else, Vmin = 0.90*Un; Vmax = 1.06*Un; end

    % Helper to compute Pmax (MW) for a single load current I (A)
    getP = @(I) bsearchP(Un,Xth,n(k),I,pfL,phi,Vmin,Vmax);

    Pmin(k) = getP(Ilo(k));
    Pmax(k) = getP(Ihi(k));
    Pcons(k)= min(Pmin(k),Pmax(k));
end

results = table(string(names)', U(:)/1e3, Pmin, Pmax, Pcons, ...
    'VariableNames',{'Point','U_kV','Pmax_at_Imin_MW','Pmax_at_Imax_MW','Pmax_conservative_MW'});
disp(results);

% ---------- helpers ----------
function Pmw = bsearchP(Un,Xth,nlines,I,pfL,phi,Vmin,Vmax)
    Sload = nlines*sqrt(3)*Un*I;               % VA
    Pload = pfL*Sload; Qload = Pload*tan(phi); % W, var (inductive)
    ok = @(Pmw) feasible(Pmw*1e6,Un,Xth,Pload,Qload,Vmin,Vmax);
    lo=0; hi=10; while ok(hi), hi=hi*2; if hi>2e4, break; end, end
    while (hi-lo)>0.1
        mid=0.5*(hi+lo); if ok(mid), lo=mid; else, hi=mid; end
    end
    Pmw=lo;
end

function tf = feasible(Pw,Un,Xth,Pload,Qload,Vmin,Vmax)
    Snet = (Pload - Pw) + 1j*(Qload - 0);      % WF pf=1 → Q=0
    V = Un;                                    % fixed-point solve
    for i=1:40
        I = conj(Snet)/(sqrt(3)*V);
        V = Un - 1j*Xth*I;
    end
    Vm = abs(V);
    tf = (Vm>=Vmin && Vm<=Vmax);
end
