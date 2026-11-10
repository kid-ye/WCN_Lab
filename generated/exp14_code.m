% Experiment 14 : Mobile Ad Hoc Network with an Omnidirectional Antenna Model
% Syllabus experiment 3.
clc; clear; close all;
rand('state',14); randn('state',14);

%% Radio parameters, omnidirectional antennas
Pt   = 20;                  % transmit power (dBm)
Gt   = 0; Gr = 0;           % omnidirectional antennas (dBi)
fc   = 2.4e9;
Sens = -85;                 % receiver sensitivity (dBm)
n_pl = 3.0;                 % path loss exponent
d0   = 1;                   % reference distance (m)
c    = 3e8;
lambda = c/fc;

PL0  = 20*log10(4*pi*d0/lambda);
Rmax = d0*10^((Pt+Gt+Gr-Sens-PL0)/(10*n_pl));
fprintf('Omnidirectional antenna model\n');
fprintf('  Transmit power        = %d dBm\n',Pt);
fprintf('  Receiver sensitivity  = %d dBm\n',Sens);
fprintf('  Path loss exponent    = %.1f\n',n_pl);
fprintf('  Reference loss at 1 m = %.2f dB\n',PL0);
fprintf('  Transmission range    = %.1f m\n\n',Rmax);

%% Network area and node counts
A      = 500;                     % square side (m)
nvec   = 10:5:60;
Ntrial = 200;

conn = zeros(size(nvec)); hops = conn; pdr = conn; deg = conn;

for k = 1:numel(nvec)
    N = nvec(k);
    cc = 0; hh = 0; pp = 0; dd = 0; cnt = 0;
    for t = 1:Ntrial
        x = A*rand(1,N); y = A*rand(1,N);
        D = sqrt((x-x').^2 + (y-y').^2);
        Adj = (D <= Rmax) & ~eye(N);
        dd = dd + mean(sum(Adj,2));

        % breadth first search from node 1 gives hop counts
        hopmat = bfs_hops(Adj);
        reach  = hopmat(1,:);
        reach(1) = 0;
        cc = cc + mean(isfinite(reach));
        fin = reach(isfinite(reach) & reach > 0);
        if ~isempty(fin), hh = hh + mean(fin); cnt = cnt + 1; end

        % connectivity based delivery ratio, 50 random source-destination
        % pairs. A pair counts as delivered when a multihop path exists; the
        % simulation does not model medium access or collisions.
        s = randi(N,1,50); d = randi(N,1,50);
        ok = 0;
        for q = 1:50
            if s(q) ~= d(q) && isfinite(hopmat(s(q),d(q))), ok = ok + 1; end
        end
        pp = pp + ok/50;
    end
    conn(k) = cc/Ntrial; hops(k) = hh/max(cnt,1);
    pdr(k)  = pp/Ntrial; deg(k)  = dd/Ntrial;
end

fprintf('%8s %14s %16s %14s %14s\n',...
        'Nodes','Mean degree','Reachable frac','Mean hops','Conn. delivery');
for k = 1:numel(nvec)
    fprintf('%8d %14.2f %16.4f %14.2f %14.4f\n',...
            nvec(k),deg(k),conn(k),hops(k),pdr(k));
end

kcrit = find(pdr > 0.95,1);
if isempty(kcrit), kcrit = numel(nvec); end
fprintf('\nThe delivery ratio first exceeds 0.95 at %d nodes, where the mean\n',nvec(kcrit));
fprintf('node degree is %.1f. Below this density the graph partitions\n',deg(kcrit));
fprintf('and packets addressed to unreachable nodes are lost, however good\n');
fprintf('the routing protocol may be.\n\n');

%% Effect of transmission range at fixed node count
N = 30;
rvec = 60:20:260;
pdr_r = zeros(size(rvec));
for k = 1:numel(rvec)
    pp = 0;
    for t = 1:Ntrial
        x = A*rand(1,N); y = A*rand(1,N);
        D = sqrt((x-x').^2 + (y-y').^2);
        Adj = (D <= rvec(k)) & ~eye(N);
        hm = bfs_hops(Adj);
        s = randi(N,1,50); d = randi(N,1,50);
        ok = 0;
        for q = 1:50
            if s(q) ~= d(q) && isfinite(hm(s(q),d(q))), ok = ok + 1; end
        end
        pp = pp + ok/50;
    end
    pdr_r(k) = pp/Ntrial;
end
fprintf('%16s %24s\n','Range (m)','Connectivity delivery ratio');
for k = 1:numel(rvec)
    fprintf('%16d %24.4f\n',rvec(k),pdr_r(k));
end
fprintf('\nRaising the transmission power extends the range and improves\n');
fprintf('connectivity, but increases interference and drains the battery,\n');
fprintf('so the useful operating point lies just above the connectivity\n');
fprintf('threshold rather than at maximum power.\n\n');

%% Figure 1 : one network realisation with routes
N = 30;
x = A*rand(1,N); y = A*rand(1,N);
D = sqrt((x-x').^2 + (y-y').^2);
Adj = (D <= Rmax) & ~eye(N);
figure; hold on;
for i = 1:N
    for j = i+1:N
        if Adj(i,j), plot([x(i) x(j)],[y(i) y(j)],'Color',[0.75 0.8 0.85]); end
    end
end
plot(x,y,'o','MarkerFaceColor',[0.20 0.40 0.65],'MarkerSize',7);
[pathN,~] = bfs_path(Adj,1,N);
if ~isempty(pathN)
    plot(x(pathN),y(pathN),'r-','LineWidth',2.5);
    plot(x(pathN),y(pathN),'rs','MarkerSize',9,'LineWidth',1.5);
end
axis equal; axis([0 A 0 A]); grid on;
xlabel('x (m)'); ylabel('y (m)');
title(sprintf('MANET, %d nodes, range %.0f m, route shown in red',N,Rmax));

%% Figure 2 : connectivity against node count
figure;
plot(nvec,pdr,'o-','LineWidth',1.8); hold on;
plot(nvec,conn,'s-','LineWidth',1.8);
grid on; axis([min(nvec) max(nvec) 0 1.05]);
xlabel('Number of nodes'); ylabel('Fraction');
legend('Connectivity delivery ratio','Reachable node fraction','Location','southeast');
title('Connectivity against Node Density');

%% Figure 3 : mean degree and hop count
figure;
subplot(1,2,1);
plot(nvec,deg,'o-','LineWidth',1.8); grid on;
xlabel('Number of nodes'); ylabel('Mean node degree');
title('Node Degree');
subplot(1,2,2);
plot(nvec,hops,'s-','LineWidth',1.8); grid on;
xlabel('Number of nodes'); ylabel('Mean hop count');
title('Route Length');

%% Figure 4 : delivery ratio against transmission range
figure;
plot(rvec,pdr_r,'o-','LineWidth',1.8); hold on;
plot([Rmax Rmax],[0 1.05],'k--','LineWidth',1.5);
grid on; axis([min(rvec) max(rvec) 0 1.05]);
xlabel('Transmission range (m)'); ylabel('Connectivity delivery ratio');
legend('Connectivity delivery ratio','Range at 20 dBm','Location','southeast');
title('Effect of Transmission Range, 30 Nodes');

function Hm = bfs_hops(Adj)
% Hop count between every pair of nodes by breadth first search.
% Unreachable pairs are returned as Inf.
N  = size(Adj,1);
Hm = inf(N);
for s = 1:N
    d = inf(1,N); d(s) = 0; q = s;
    while ~isempty(q)
        v = q(1); q(1) = [];
        nb = find(Adj(v,:) & isinf(d));
        d(nb) = d(v) + 1;
        q = [q nb];
    end
    Hm(s,:) = d;
end
end

function [p,h] = bfs_path(Adj,s,t)
% Shortest hop path from s to t by breadth first search.
N = size(Adj,1);
prev = zeros(1,N); d = inf(1,N); d(s) = 0; q = s;
while ~isempty(q)
    v = q(1); q(1) = [];
    nb = find(Adj(v,:) & isinf(d));
    d(nb) = d(v) + 1; prev(nb) = v;
    q = [q nb];
end
if isinf(d(t)), p = []; h = Inf; return; end
p = t; while p(1) ~= s, p = [prev(p(1)) p]; end
h = d(t);
end
