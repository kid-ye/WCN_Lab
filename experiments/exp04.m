% Experiment 4 : Log-Normal Shadowing and Rayleigh Fading
clc; clear; close all;
randn('state',4); rand('state',4);

%% Parameters
fc     = 900e6;         % carrier frequency (Hz)
v      = 60;            % vehicle speed (km/h)
sigma  = 8;             % shadowing standard deviation (dB)
M      = 32;            % number of oscillators
c      = 3e8;
EIRP   = 55;  Sens = -104;  maxPL = EIRP - Sens;
L1_sub = 116.46;        % suburban Hata path loss at 1 km, from Experiment 3
nEff   = 3.5225;

fm = (v/3.6)*fc/c;
fprintf('Maximum Doppler shift = %.2f Hz\n\n',fm);

%% Part A : log-normal shadowing
X = sigma*randn(1,2e6);
fprintf('Shadowing standard deviation (simulated) = %.4f dB\n',std(X));

target = [0.75 0.90 0.95 0.99];
fprintf('\n%16s %14s %16s\n','Edge coverage','Margin (dB)','Simulated');
for k = 1:numel(target)
    Mk = sigma*sqrt(2)*erfcinv(2*(1-target(k)));
    fprintf('%15.0f%% %14.4f %15.4f\n',100*target(k),Mk,mean(X > -Mk));
end

M90 = sigma*sqrt(2)*erfcinv(2*0.10);
R0  = 10^((maxPL       - L1_sub)/(10*nEff));
R1  = 10^((maxPL - M90 - L1_sub)/(10*nEff));
fprintf('\nCoverage radius on median path loss = %.2f km\n',R0);
fprintf('Coverage radius with %.2f dB margin = %.2f km\n',M90,R1);
fprintf('Radius reduced by a factor of %.3f, sites increased by %.3f\n\n',R0/R1,(R0/R1)^2);

%% Part B : Rayleigh fading by the sum of sinusoids
fs   = 50*fm;                    % sampling rate
Tper = 2000;                     % record length in Doppler periods
N    = round(Tper*fs/fm);
t    = (0:N-1)'/fs;
n    = 1:M;

theta = 2*pi*rand - pi;
phi   = 2*pi*rand(1,M) - pi;
psi   = 2*pi*rand(1,M) - pi;
alpha = (2*pi*n - pi + theta)/(4*M);

arg = 2*pi*fm*(t*cos(alpha)) + repmat(phi,N,1);
hc  = sqrt(2/M)*(cos(arg)*cos(psi)');
hs  = sqrt(2/M)*(cos(arg)*sin(psi)');
h   = hc + 1j*hs;
r   = abs(h);
rms = sqrt(mean(r.^2));

fprintf('Mean power E|h|^2 = %.4f  (theory 1.0000)\n',mean(r.^2));
fprintf('Envelope mean     = %.4f  (theory %.4f)\n',mean(r),sqrt(pi)/2);
fprintf('E[r^4]/E[r^2]^2   = %.4f  (theory 2.0000)\n\n',mean(r.^4)/mean(r.^2)^2);

%% Autocorrelation against Clarke theory
maxlag = round(3*fs/fm);
hz  = h - mean(h);
Rxx = zeros(maxlag+1,1);
for k = 0:maxlag
    Rxx(k+1) = mean(conj(hz(1:end-k)).*hz(1+k:end));
end
Rxx = real(Rxx)/real(Rxx(1));
tau = (0:maxlag)'/fs;
J0  = besselj(0,2*pi*fm*tau);
fprintf('Maximum |ACF - J0(2*pi*fm*tau)| = %.4f\n\n',max(abs(Rxx-J0)));

%% Level crossing rate and average fade duration
rho_dB = [-20 -15 -10 -5 0 3];
fprintf('%8s %12s %12s %14s %14s\n','rho(dB)','NR sim','NR theory','AFD sim (ms)','AFD thy (ms)');
NRs = zeros(size(rho_dB)); AFDs = NRs; NRt = NRs; AFDt = NRs;
for k = 1:numel(rho_dB)
    rho   = 10^(rho_dB(k)/20);
    level = rho*rms;
    below = r < level;
    up    = sum(~below(2:end) & below(1:end-1));
    NRs(k)  = up/(N/fs);
    AFDs(k) = (sum(below)/fs)/up;
    NRt(k)  = sqrt(2*pi)*fm*rho*exp(-rho^2);
    AFDt(k) = (exp(rho^2)-1)/(rho*fm*sqrt(2*pi));
    fprintf('%8.0f %12.2f %12.2f %14.3f %14.3f\n',...
            rho_dB(k),NRs(k),NRt(k),AFDs(k)*1e3,AFDt(k)*1e3);
end

fprintf('\nGSM frame duration = 4.615 ms\n');
for vp = [5 60]
    fmp = (vp/3.6)*fc/c;
    rho = 10^(-10/20);
    afd = (exp(rho^2)-1)/(rho*fmp*sqrt(2*pi));
    fprintf('At %2d km/h : fm = %5.2f Hz, AFD at -10 dB = %7.3f ms\n',vp,fmp,afd*1e3);
end

%% Figure 1 : fading envelope
ns = round(20*fs/fm);
figure;
plot(t(1:ns)*1000,20*log10(r(1:ns)/rms),'LineWidth',1); hold on;
plot([0 t(ns)*1000],[-10 -10],'r--','LineWidth',1.5);
grid on; xlabel('Time (ms)'); ylabel('Envelope relative to rms (dB)');
legend('Rayleigh envelope','-10 dB threshold','Location','southwest');
title(sprintf('Fading Envelope at %d km/h, f_m = %.0f Hz',v,fm));

%% Figure 2 : envelope distribution
figure;
[cnt,ctr] = hist(r/rms,80);
bar(ctr,cnt/(sum(cnt)*(ctr(2)-ctr(1))),1); hold on;
xr = linspace(0,max(ctr),300);
plot(xr,2*xr.*exp(-xr.^2),'r','LineWidth',2);
grid on; xlabel('Envelope / rms'); ylabel('Probability Density');
legend('Simulated','Rayleigh pdf'); title('Envelope Distribution');

%% Figure 3 : Doppler spectrum
seg = 4096; nseg = floor(N/seg);
w = 0.5 - 0.5*cos(2*pi*(0:seg-1)'/(seg-1));
P = zeros(seg,1);
for k = 1:nseg
    P = P + abs(fft(h((k-1)*seg+1:k*seg).*w)).^2;
end
P = fftshift(P/nseg);
fax = (-seg/2:seg/2-1)'*fs/seg;
figure;
plot(fax,10*log10(P/max(P)),'LineWidth',1.2); hold on;
plot([-fm -fm],[-60 5],'r--',[fm fm],[-60 5],'r--','LineWidth',1.5);
grid on; axis([-2*fm 2*fm -60 5]);
xlabel('Frequency (Hz)'); ylabel('Normalised PSD (dB)');
title('Doppler Spectrum');

%% Figure 4 : level crossing rate and fade duration
figure;
subplot(1,2,1);
semilogy(rho_dB,NRs,'o-','LineWidth',1.6); hold on;
semilogy(rho_dB,NRt,'k--','LineWidth',1.6);
grid on; xlabel('\rho (dB)'); ylabel('Crossings per second');
legend('Simulated','Theory','Location','southeast'); title('Level Crossing Rate');
subplot(1,2,2);
semilogy(rho_dB,AFDs*1e3,'o-','LineWidth',1.6); hold on;
semilogy(rho_dB,AFDt*1e3,'k--','LineWidth',1.6);
grid on; xlabel('\rho (dB)'); ylabel('Duration (ms)');
legend('Simulated','Theory','Location','southeast'); title('Average Fade Duration');
