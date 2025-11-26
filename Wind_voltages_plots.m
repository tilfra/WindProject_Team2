clear; clc; close all

% Data for points
names={'10 kV point 1','10 kV bus','20/10 kV side (20 kV)','20 kV bus'};
U_kV=[10 10 20 20];
Sk_MVA=[500 900 1100 1300];
Imin=[150 150 130 130];
Imax=[370 390 345 345];
n=[2 4 2 2];

% 10 kV
idx = [1 2];

pfL = 0.9; phi = acos(pfL);
RX  = 0; % R/X ratio
Pmw = 0:1:100; % sweep

for p = 1:numel(idx)
    k=idx(p);
    Un=U_kV(k)*1e3;
    Sk=Sk_MVA(k)*1e6;
    Zm=(Un^2)/Sk;
    X=Zm/sqrt(1+RX^2);
    %X=0.3;
    R=RX*X;
    Z=R+1j*X;

    % Voltage limits
    Vmin = 0.90*Un;  
    Vmax = 1.06*Un;

    % Load powers for min/max load
    Smin = n(k)*sqrt(3)*Un*Imin(k);   Pmin = pfL*Smin; Qmin = Pmin*tan(phi);
    Smax = n(k)*sqrt(3)*Un*Imax(k);   Pmax = pfL*Smax; Qmax = Pmax*tan(phi);

    Vm_min = zeros(size(Pmw));
    Vm_max = zeros(size(Pmw));

    for i = 1:numel(Pmw)
        Pw = Pmw(i)*1e6;

        % Min load curve
        V = Un; 
        for it = 1:40
            Iline = conj((Pmin - Pw) + 1j*Qmin) / (sqrt(3)*V);
            V = Un - Z*Iline;
        end
        Vm_min(i) = abs(V)/1e3;

        % Max load curve
        V = Un;
        for it = 1:40
            Iline = conj((Pmax - Pw) + 1j*Qmax) / (sqrt(3)*V);
            V = Un - Z*Iline;
        end
        Vm_max(i) = abs(V)/1e3;
    end

    % Plot
    figure('Color','w'); hold on; grid on; box on
    yline(Vmin/1e3,'k.-.','DisplayName','Min voltage')
    yline(Vmax/1e3,'k.-','DisplayName','Max voltage')
    plot(Pmw, Vm_min, 'LineWidth',1.6, 'DisplayName','Min load')
    plot(Pmw, Vm_max, 'LineWidth',1.6, 'DisplayName','Max load')
    xlabel('Wind power P_{wind} (MW)')
    ylabel('|U_{PCC}| (kV)')
    title(sprintf('%s  (U=%.0f kV, S_k=%.0f MVA, R/X=%.f)',names{k}, Un/1e3, Sk/1e6, RX));
    legend('Location','southwest'); xlim([0 100]);
end
