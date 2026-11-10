% Experiment 1 : Frequency Reuse, Co-channel Interference and Trunked Capacity
clc; clear; close all;

%% Parameters
BW_total = 25e6;          % total duplex spectrum (Hz)
BW_ch    = 200e3;         % channel spacing (Hz)
n        = 4;             % path loss exponent
i0       = 6;             % first-tier co-channel interferers
SI_req   = 18;            % required S/I (dB)
GoS      = 0.02;          % grade of service
E_sub    = 0.025;         % traffic per subscriber (erlang)
Area     = 400;           % service area (km^2)
Density  = 800;           % subscribers per km^2

Ctot = floor(BW_total/BW_ch);
Nset = [3 4 7 9 12];

%% Reuse ratio, S/I, channels per cell, capacity
R = zeros(numel(Nset),7);
for k = 1:numel(Nset)
    N  = Nset(k);
    DR = sqrt(3*N);                       % co-channel reuse ratio
    SI = 10*log10(DR^n / i0);             % signal to interference ratio (dB)
    C  = floor(Ctot/N);                   % channels per cell
    A  = erlangb_traffic(C,GoS);          % offered traffic at 2 % blocking
    Ac = A*(1-GoS);                       % carried traffic
    S  = A/E_sub;                         % subscribers per cell
    R(k,:) = [N DR SI C A Ac S];
end

fprintf('Total duplex channels = %d\n\n',Ctot);
fprintf('  N     D/R    S/I(dB)   Ch/cell   Offered(E)   Carried(E)   Subs/cell\n');
fprintf(' %2d  %6.3f  %7.2f  %7d  %10.3f  %10.3f  %9.0f\n',R.');

feas = R(R(:,3)>=SI_req,:);
Nmin = feas(1,1);  Cmin = feas(1,4);  Smin = feas(1,7);
fprintf('\nMinimum feasible cluster size N = %d\n',Nmin);
fprintf('Channels per cell               = %d\n',Cmin);
fprintf('Subscribers per cell            = %.0f\n',Smin);

%% Trunking efficiency
i7 = find(R(:,1)==7); i12 = find(R(:,1)==12);
fprintf('\nN = 7 to N = 12 : channels fall %.1f %%, offered traffic falls %.1f %%\n',...
        100*(1-R(i12,4)/R(i7,4)), 100*(1-R(i12,5)/R(i7,5)));

%% Cell count for the given subscriber density
Ncell = ceil(Density*Area/Smin);
Acell = Area/Ncell;
Rcell = sqrt(2*Acell/(3*sqrt(3)));
fprintf('\nSubscribers to be served = %d\n',Density*Area);
fprintf('Cells required           = %d\n',Ncell);
fprintf('Cell radius              = %.3f km\n',Rcell);

%% Sectoring
fprintf('\nSectoring:\n');
for ii = [6 2 1]
    for N = [1 3 4 7 9 12 13 16 19]
        if 10*log10((3*N)^(n/2)/ii) >= SI_req, break; end
    end
    fprintf('  i0 = %d : N = %2d, channels/cell = %d\n',ii,N,floor(Ctot/N));
end

%% Figure 1 : frequency group assignment for N = 7
Rc = 1; q = []; r = [];
for a = -3:3
    for b = -3:3
        if abs(a+b) <= 3, q(end+1,1) = a; r(end+1,1) = b; end
    end
end
xc = sqrt(3)*Rc*(q + r/2);  yc = 1.5*Rc*r;
grp = mod(q + 3*r,7) + 1;
th  = (0:5)'*pi/3 + pi/6;
col = [0.85 0.33 0.31;0.33 0.55 0.78;0.47 0.72 0.42;0.95 0.75 0.30;...
       0.62 0.45 0.75;0.45 0.75 0.78;0.80 0.55 0.40];
figure; hold on;
for k = 1:numel(xc)
    patch(xc(k)+Rc*cos(th), yc(k)+Rc*sin(th), col(grp(k),:));
    text(xc(k),yc(k),sprintf('%d',grp(k)),'HorizontalAlignment','center');
end
axis equal off; title('Frequency Group Assignment, N = 7');

%% Figure 2 : S/I and capacity against cluster size
figure;
subplot(1,2,1);
plot(R(:,1),R(:,3),'o-','LineWidth',2); hold on;
plot([min(Nset) max(Nset)],[SI_req SI_req],'k--','LineWidth',1.5);
grid on; xlabel('Cluster size N'); ylabel('S/I (dB)');
legend('S/I','Requirement','Location','southeast');
title('Interference Performance');
subplot(1,2,2);
plot(R(:,1),R(:,7),'s-','LineWidth',2);
grid on; xlabel('Cluster size N'); ylabel('Subscribers per cell');
title('Traffic Capacity');

function Pb = erlangb(C,A)
% Blocking probability of a lost calls cleared system, stable recursion.
Pb = 1;
for k = 1:C
    Pb = (A*Pb)/(k + A*Pb);
end
end

function A = erlangb_traffic(C,Pt)
% Offered traffic carried by C channels at blocking probability Pt.
if C == 0, A = 0; return; end
lo = 0; hi = 2*C;
while erlangb(C,hi) < Pt, hi = 2*hi; end
for k = 1:200
    mid = 0.5*(lo+hi);
    if erlangb(C,mid) < Pt, lo = mid; else hi = mid; end
end
A = 0.5*(lo+hi);
end
