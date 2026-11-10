% Experiment 13 : Wireless LAN Basic Service Set - CSMA/CA and Scheduled Access
% Syllabus experiments 1 and 2 : packet transmission performance of a BSS,
% and comparison of Wi-Fi with WiMAX.
clc; clear; close all;
rand('state',13);

%% IEEE 802.11 DCF parameters
Rdata   = 1e6;            % channel bit rate (bit/s)
Payload = 8184;           % payload (bits)
MAChdr  = 272;
PHYhdr  = 128;
ACK     = 112 + PHYhdr;
slot    = 20e-6;
SIFS    = 28e-6;
DIFS    = 128e-6;
delta   = 1e-6;           % propagation delay
W       = 32;             % minimum contention window
m       = 3;              % maximum backoff stage

H  = MAChdr + PHYhdr;
Ts = (H + Payload)/Rdata + SIFS + delta + ACK/Rdata + DIFS + delta;
Tc = (H + Payload)/Rdata + DIFS + delta;

fprintf('IEEE 802.11 DCF parameters\n');
fprintf('  Channel bit rate              = %.1f Mbit/s\n',Rdata/1e6);
fprintf('  Payload                       = %d bits\n',Payload);
fprintf('  Successful transmission Ts= %.2f ms\n',Ts*1e3);
fprintf('  Collision duration Tc         = %.2f ms\n',Tc*1e3);
fprintf('  Minimum contention window = %d\n',W);
fprintf('  Idle slot                     = %.0f us\n\n',slot*1e6);

%% Saturation throughput, Bianchi analysis
nvec = [2 5 10 15 20 30 40 50];
S_an = zeros(size(nvec));
tau_v = zeros(size(nvec));
for k = 1:numel(nvec)
    n = nvec(k);
    p = 0.1;
    for it = 1:500
        tau = 2*(1-2*p)/((1-2*p)*(W+1) + p*W*(1-(2*p)^m));
        pn  = 1 - (1-tau)^(n-1);
        p   = p + 0.5*(pn - p);
    end
    tau = 2*(1-2*p)/((1-2*p)*(W+1) + p*W*(1-(2*p)^m));
    Ptr = 1 - (1-tau)^n;
    Ps  = n*tau*(1-tau)^(n-1)/Ptr;
    S_an(k) = Ps*Ptr*Payload/Rdata / ...
              ((1-Ptr)*slot + Ptr*Ps*Ts + Ptr*(1-Ps)*Tc);
    tau_v(k) = tau;
end

%% Saturation throughput, slotted simulation
Nslots = 2e5;
S_sim  = zeros(size(nvec));
Pcoll  = zeros(size(nvec));
for k = 1:numel(nvec)
    n  = nvec(k);
    bo = floor(rand(1,n)*W);
    st = zeros(1,n);
    T  = 0; succ = 0; coll = 0;
    for it = 1:Nslots
        ready = (bo == 0);
        nr = sum(ready);
        if nr == 0
            T = T + slot; bo = bo - 1;
        elseif nr == 1
            T = T + Ts; succ = succ + 1;
            st(ready) = 0; bo(ready) = floor(rand(1,nr)*W);
            bo(~ready) = bo(~ready) - 1;
        else
            T = T + Tc; coll = coll + 1;
            st(ready) = min(st(ready)+1,m);
            bo(ready) = floor(rand(1,nr).*(W*2.^st(ready)));
            bo(~ready) = bo(~ready) - 1;
        end
        bo = max(bo,0);
    end
    S_sim(k) = succ*Payload/Rdata/T;
    Pcoll(k) = coll/(succ+coll);
end

%% WiMAX scheduled access, no contention
Wpay  = 8184;
Wover = 48 + 88;                         % MAC header and preamble overhead (bits)
Wrate = 1e6;
S_wmx = Wpay/(Wpay+Wover)*ones(size(nvec));

fprintf('%8s %14s %14s %14s %14s\n',...
        'Stations','DCF analysis','DCF simulation','Collision prob','WiMAX TDMA');
for k = 1:numel(nvec)
    fprintf('%8d %14.4f %14.4f %14.4f %14.4f\n',...
            nvec(k),S_an(k),S_sim(k),Pcoll(k),S_wmx(k));
end

fprintf('\nThe normalised throughput of the contention based DCF falls as\n');
fprintf('stations are added, because the collision probability rises. The\n');
fprintf('scheduled WiMAX access retains a constant efficiency because the\n');
fprintf('base station allocates the slots and no collisions occur. The\n');
fprintf('penalty is the signalling required to request and grant slots.\n\n');

%% Packet delivery ratio against offered load, non-saturated
lam = 0.05:0.05:1.2;
n   = 10;
pdr = zeros(size(lam));
for k = 1:numel(lam)
    G   = lam(k);
    Ps  = exp(-2*G);                            % slotted access approximation
    pdr(k) = min(1,Ps + (1-Ps)*0.5);
end

%% Figure 1 : throughput against number of stations
figure;
plot(nvec,S_an,'o-','LineWidth',1.8); hold on;
plot(nvec,S_sim,'s--','LineWidth',1.8);
plot(nvec,S_wmx,'d-','LineWidth',1.8);
grid on; axis([0 50 0 1]);
xlabel('Number of stations'); ylabel('Normalised throughput');
legend('802.11 DCF, analysis','802.11 DCF, simulation','WiMAX scheduled access',...
       'Location','east');
title('Basic Service Set Throughput against Load');

%% Figure 2 : collision probability
figure;
plot(nvec,Pcoll,'o-','LineWidth',1.8); hold on;
plot(nvec,tau_v,'s-','LineWidth',1.8);
grid on; xlabel('Number of stations'); ylabel('Probability');
legend('Collision probability','Transmission probability \tau','Location','east');
title('Contention Behaviour of the Distributed Coordination Function');

%% Figure 3 : time spent in each state
figure;
Tsucc = zeros(size(nvec)); Tcol = Tsucc; Tidl = Tsucc;
for k = 1:numel(nvec)
    n = nvec(k); tau = tau_v(k);
    Ptr = 1-(1-tau)^n; Ps = n*tau*(1-tau)^(n-1)/Ptr;
    tot = (1-Ptr)*slot + Ptr*Ps*Ts + Ptr*(1-Ps)*Tc;
    Tidl(k)  = (1-Ptr)*slot/tot;
    Tsucc(k) = Ptr*Ps*Ts/tot;
    Tcol(k)  = Ptr*(1-Ps)*Tc/tot;
end
bar(nvec,[Tsucc' Tcol' Tidl'],'stacked');
xlabel('Number of stations'); ylabel('Fraction of channel time');
legend('Successful transmission','Collision','Idle','Location','east');
title('Channel Time Utilisation');
