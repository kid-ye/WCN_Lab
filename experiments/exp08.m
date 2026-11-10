% Experiment 8 : GMSK Modulation and Spectral Analysis
clc; clear; close all;
rand('state',8);

%% Parameters
Rb   = 270833.33;     % GSM gross bit rate (bit/s)
sps  = 16;            % samples per bit
BT   = 0.3;           % Gaussian bandwidth-time product
span = 4;             % Gaussian pulse span (bit periods)
fs   = sps*Rb;
Tb   = 1/Rb;

fprintf('Bit rate          = %.2f kbit/s\n',Rb/1000);
fprintf('Bit period        = %.4f us\n',Tb*1e6);
fprintf('Sampling rate     = %.4f MHz\n',fs/1e6);
fprintf('BT product        = %.2f\n\n',BT);

%% Gaussian frequency pulse
t     = (-span:1/sps:span)';
sigma = sqrt(log(2))/(2*pi*BT);
gauss = exp(-t.^2/(2*sigma^2))/(sqrt(2*pi)*sigma);
rect  = ones(sps,1)/sps;
gp    = conv(gauss,rect,'same');
gp    = gp/sum(gp);                      % unit area, each bit adds pi/2

fprintf('Gaussian pulse: sigma = %.4f bit periods, area = %.6f\n\n',sigma,sum(gp));

%% GSM normal burst
tsc   = [0 0 1 0 0 1 0 1 1 1 0 0 0 0 1 0 0 0 1 0 0 1 0 1 1 1];
burst = [0 0 0 randi([0 1],1,57) 0 tsc 0 randi([0 1],1,57) 0 0 0];
fprintf('Normal burst length = %d bits\n',numel(burst));
fprintf('Burst duration      = %.4f ms\n',(numel(burst)+8.25)*Tb*1e3);
fprintf('TDMA frame duration = %.4f ms\n\n',8*(numel(burst)+8.25)*Tb*1e3);

%% Modulate the burst
a  = 2*burst - 1;
up = zeros(1,numel(a)*sps); up(1:sps:end) = a;
fq = conv(up,gp','same');
ph = (pi/2)*cumsum(fq);
s  = exp(1j*ph);

fprintf('Envelope variation max|s| - min|s| = %.3e\n',max(abs(s))-min(abs(s)));
fprintf('GMSK is a constant envelope modulation.\n\n');

%% Spectra for several BT products
Nb   = 4000;
bits = randi([0 1],1,Nb);
BTs  = [0.2 0.3 0.5 1e3];
lab  = {'BT = 0.2','BT = 0.3 (GSM)','BT = 0.5','MSK, no filter'};
nfft = 2048;
P    = cell(1,4);
fprintf('%-18s %14s %12s\n','Filter','99 % BW (kHz)','BW / Rb');
for k = 1:4
    sg = sqrt(log(2))/(2*pi*BTs(k));
    gg = exp(-t.^2/(2*sg^2))/(sqrt(2*pi)*sg);
    gk = conv(gg,rect,'same'); gk = gk/sum(gk);
    ak = 2*bits - 1;
    uk = zeros(1,numel(ak)*sps); uk(1:sps:end) = ak;
    sk = exp(1j*(pi/2)*cumsum(conv(uk,gk','same')));

    w = 0.5 - 0.5*cos(2*pi*(0:nfft-1)'/(nfft-1));
    nseg = floor(numel(sk)/nfft); Pk = zeros(nfft,1);
    for m = 1:nseg
        Pk = Pk + abs(fft(sk((m-1)*nfft+1:m*nfft)'.*w)).^2;
    end
    P{k} = fftshift(Pk/nseg);
    fax  = (-nfft/2:nfft/2-1)'*fs/nfft;

    tot = sum(P{k}); ctr = nfft/2+1; acc = P{k}(ctr); j = 0;
    while acc < 0.99*tot && ctr-j-1 >= 1 && ctr+j+1 <= nfft
        j = j + 1; acc = acc + P{k}(ctr-j) + P{k}(ctr+j);
    end
    bw = 2*j*(fax(2)-fax(1));
    fprintf('%-18s %14.1f %12.2f\n',lab{k},bw/1000,bw/Rb);
end
fprintf('\nThe GSM channel spacing is 200 kHz. The occupied bandwidth exceeds\n');
fprintf('this because channel spacing is set by the tolerable adjacent\n');
fprintf('channel power, not by the occupied bandwidth.\n');

%% Figure 1 : frequency pulse
figure;
plot((0:numel(gp)-1)/sps - span,gp*sps,'LineWidth',2); hold on;
plot([-0.5 -0.5 0.5 0.5],[0 1 1 0],'k--','LineWidth',1.5);
grid on; xlabel('Time (bit periods)'); ylabel('Normalised Amplitude');
legend('Gaussian pulse, BT = 0.3','Rectangular pulse, MSK');
title('Frequency Pulse Shapes');

%% Figure 2 : complex envelope
tt = (0:numel(s)-1)/fs;
figure;
plot(tt*1e6,real(s),'LineWidth',1.2); hold on;
plot(tt*1e6,imag(s),'LineWidth',1.2);
plot(tt*1e6,abs(s),'k','LineWidth',1.8);
grid on; axis([0 60 -1.3 1.3]);
xlabel('Time (\mus)'); ylabel('Amplitude');
legend('In-phase','Quadrature','Envelope','Location','southeast');
title('GMSK Complex Envelope');

%% Figure 3 : spectra
figure; hold on;
for k = 1:4
    plot(fax/1000,10*log10(P{k}/max(P{k})),'LineWidth',1.6);
end
plot([100 100],[-90 5],'k--',[-100 -100],[-90 5],'k--','LineWidth',1.4);
grid on; axis([-400 400 -90 5]);
xlabel('Frequency Offset (kHz)'); ylabel('Normalised PSD (dB)');
legend([lab,{'200 kHz channel edge'}]);
title('Effect of the BT Product on the Transmitted Spectrum');

%% Figure 4 : phase trajectory
figure;
plot((0:numel(s)-1)/sps,unwrap(angle(s))/pi,'LineWidth',1.6);
grid on; xlabel('Time (bit periods)'); ylabel('Phase / \pi');
title('Phase Trajectory : Each Bit Advances the Phase by \pm\pi/2');
