% Experiment 12 : 5G NR Numerology and Waveform Comparison
clc; clear; close all;
rand('state',12); randn('state',12);

%% Numerology definition
mu   = 0:4;
scs  = 15e3*2.^mu;                    % subcarrier spacing (Hz)
Nfft = 4096;

fprintf('5G NR numerology, subcarrier spacing = 15 x 2^mu kHz\n\n');
fprintf('%5s %10s %12s %12s %12s %12s %10s\n',...
        'mu','SCS (kHz)','Tu (us)','TCP (us)','Tsym (us)','Slot (ms)','Slots/ms');
Tu = zeros(size(mu)); Tcp = Tu; Tsym = Tu; Tslot = Tu;
for k = 1:numel(mu)
    Tu(k)    = 1/scs(k);
    Tcp(k)   = 144/2048*Tu(k);                 % normal cyclic prefix
    Tsym(k)  = Tu(k) + Tcp(k);
    Tslot(k) = 14*Tsym(k);
    fprintf('%5d %10.0f %12.3f %12.3f %12.3f %12.4f %10.0f\n',...
            mu(k),scs(k)/1000,Tu(k)*1e6,Tcp(k)*1e6,Tsym(k)*1e6,Tslot(k)*1e3,1e-3/Tslot(k));
end

%% Delay spread and Doppler limits
fprintf('\n%5s %10s %18s %18s %16s\n',...
        'mu','SCS (kHz)','Max delay spread','Max Doppler','Max speed at 3.5 GHz');
fc = 3.5e9; c = 3e8;
for k = 1:numel(mu)
    tmax = Tcp(k);
    % Design criterion adopted for this experiment: the Doppler shift is
    % held below one per cent of the subcarrier spacing so that the
    % intercarrier interference remains small. This is a working heuristic
    % for the exercise, not a value specified by 3GPP.
    fdmax = 0.01*scs(k);
    vmax  = fdmax*c/fc*3.6;
    fprintf('%5d %10.0f %15.3f us %15.0f Hz %13.0f km/h\n',...
            mu(k),scs(k)/1000,tmax*1e6,fdmax,vmax);
end
fprintf('\nA wider subcarrier spacing shortens the symbol and the cyclic\n');
fprintf('prefix, which reduces latency and raises the tolerable Doppler\n');
fprintf('shift, but reduces the delay spread that can be absorbed.\n\n');

%% Deployment mapping
% A numerology is not tied to a frequency range by the standard. The
% association below is a typical deployment choice, arrived at from the
% delay spread and mobility of each scenario, not a rule.
fprintf('%5s %10s %-30s %s\n','mu','SCS (kHz)','Typical deployment','Governing constraint');
dep = {'Wide area coverage, large cells','Large delay spread';
       'Urban macro and micro cells','Balanced';
       'Small cells, low latency services','Latency and Doppler';
       'Dense small cell and millimetre wave','Phase noise and Doppler';
       'Wideband millimetre wave','Phase noise and Doppler'};
for k = 1:numel(mu)
    fprintf('%5d %10.0f %-30s %s\n',mu(k),scs(k)/1000,dep{k,1},dep{k,2});
end
fprintf('\nA given numerology may be used in more than one frequency band.\n');
fprintf('The band and the numerology are linked through the delay spread\n');
fprintf('and mobility of the deployment, not by any direct rule.\n');

%% Waveform generation and spectral comparison
Nrb = 51; Nsc = 12;
fprintf('\n%5s %10s %14s %16s %16s\n',...
        'mu','SCS (kHz)','Occupied BW','Samples per slot','PAPR 99 % (dB)');
for k = 1:3
    Nu   = Nrb*Nsc;
    Nf   = 2^ceil(log2(Nu/0.85));
    used = [Nf/2-Nu/2+1:Nf/2, Nf/2+2:Nf/2+1+Nu/2];
    Ns   = 14;
    b    = randi([0 1],2,Nu*Ns);
    X    = reshape(((1-2*b(1,:))+1j*(1-2*b(2,:)))/sqrt(2),Nu,Ns);
    G    = zeros(Nf,Ns); G(used,:) = X;
    x    = ifft(ifftshift(G,1),Nf)*sqrt(Nf);
    ncp  = round(144/2048*Nf);
    xc   = [x(end-ncp+1:end,:); x];
    pp   = 10*log10(max(abs(xc).^2,[],1)./mean(abs(xc).^2,1));
    fprintf('%5d %10.0f %11.2f MHz %16d %16.2f\n',...
            mu(k),scs(k)/1000,Nu*scs(k)/1e6,numel(xc),prctile_local(pp,99));
end

%% Figure 1 : symbol duration and slot duration
figure;
subplot(1,2,1);
semilogy(mu,Tsym*1e6,'o-','LineWidth',2); grid on;
xlabel('\mu'); ylabel('Symbol duration (\mus)');
set(gca,'XTick',mu); title('OFDM Symbol Duration');
subplot(1,2,2);
semilogy(mu,Tslot*1e3,'s-','LineWidth',2); grid on;
xlabel('\mu'); ylabel('Slot duration (ms)');
set(gca,'XTick',mu); title('Slot Duration');

%% Figure 2 : cyclic prefix against delay spread
figure;
bar(mu,Tcp*1e6); hold on;
plot([-0.5 4.5],[0.5 0.5],'r--','LineWidth',2);
plot([-0.5 4.5],[0.1 0.1],'k--','LineWidth',2);
grid on; set(gca,'XTick',mu);
xlabel('\mu'); ylabel('Cyclic prefix duration (\mus)');
legend('Normal cyclic prefix','Typical urban delay spread',...
       'Typical indoor delay spread');
title('Cyclic Prefix against Channel Delay Spread');

%% Figure 3 : subcarrier spectra for three numerologies
figure; hold on;
for k = 1:3
    f = linspace(-4*scs(k),4*scs(k),2000);
    plot(f/1e3,20*log10(abs(sinc(f/scs(k)))+1e-6),'LineWidth',1.6);
end
grid on; axis([-500 500 -60 5]);
xlabel('Frequency (kHz)'); ylabel('Normalised amplitude (dB)');
legend('15 kHz','30 kHz','60 kHz');
title('Subcarrier Spectra for Three Numerologies');

%% Figure 4 : latency against Doppler tolerance
figure;
plotyy_x = mu;
[ax,h1,h2] = plotyy(plotyy_x,Tslot*1e3,plotyy_x,0.01*scs/1000);
set(h1,'Marker','o','LineWidth',2); set(h2,'Marker','s','LineWidth',2);
set(ax,'XTick',mu);
xlabel('\mu'); ylabel(ax(1),'Slot duration (ms)'); ylabel(ax(2),'Tolerable Doppler (kHz)');
legend('Slot duration','Tolerable Doppler');
title('Latency against Doppler Tolerance');
grid on;

function y = prctile_local(x,p)
x = sort(x(:));
n = numel(x);
y = interp1(100*((1:n)'-0.5)/n, x, min(max(p,100*0.5/n),100*(n-0.5)/n));
end
