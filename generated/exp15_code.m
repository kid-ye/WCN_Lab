% Experiment 15 : Optical Fibre Link - Analog Link, Digital Link and Fibre Loss
% Syllabus experiments 4, 5 and 6, implemented as MATLAB link computations.
clc; clear; close all;

%% PART A : Fibre attenuation against wavelength
lam  = 0.8:0.005:1.7;                          % wavelength (um)
Aray = 0.9./lam.^4;                            % Rayleigh scattering
Air  = 6.5e-3*exp((lam-1.55)/0.06);            % infrared absorption
Aoh  = 0.35*exp(-((lam-1.39)/0.025).^2);       % hydroxyl ion peak
Atot = Aray + Air + Aoh + 0.02;

fprintf('PART A : Fibre attenuation\n');
win = [0.85 1.31 1.55];
nam = {'First','Second','Third'};
fprintf('%10s %14s %18s\n','Window','Wavelength (um)','Attenuation (dB/km)');
alp = zeros(1,3);
for k = 1:3
    [~,i] = min(abs(lam-win(k)));
    alp(k) = Atot(i);
    fprintf('%10s %14.2f %18.4f\n',nam{k},win(k),alp(k));
end
[amin,imin] = min(Atot);
fprintf('\nMinimum attenuation %.4f dB/km at %.3f um\n',amin,lam(imin));
fprintf('Rayleigh scattering varies as 1/lambda^4 and dominates below\n');
fprintf('1.3 um; infrared absorption rises beyond 1.6 um. The minimum\n');
fprintf('between them is the reason long haul systems work at 1.55 um.\n\n');

%% PART B : Optical power budget and maximum link length
Ptx     = 0;                % launched power (dBm)
Prx     = -30;              % receiver sensitivity (dBm)
Lconn   = 0.5;              % connector loss (dB each)
Nconn   = 2;
Lsplice = 0.1;              % splice loss (dB each)
Lspan   = 2;                % splice spacing (km)
Margin  = 6;                % system margin (dB)

fprintf('PART B : Power budget at 1.55 um\n');
fprintf('  Launched power        = %.1f dBm\n',Ptx);
fprintf('  Receiver sensitivity  = %.1f dBm\n',Prx);
fprintf('  Total budget          = %.1f dB\n',Ptx-Prx);
fprintf('  Connector loss        = %.1f dB\n',Nconn*Lconn);
fprintf('  System margin         = %.1f dB\n',Margin);
avail = (Ptx-Prx) - Nconn*Lconn - Margin;
fprintf('  Available for fibre   = %.1f dB\n',avail);
Lmax = avail/(alp(3) + Lsplice/Lspan);
fprintf('  Fibre attenuation     = %.4f dB/km\n',alp(3));
fprintf('  Splice loss           = %.3f dB/km\n',Lsplice/Lspan);
fprintf('  Maximum link length   = %.2f km\n\n',Lmax);

%% PART C : Digital link, bit error rate against received power
Rb    = 2.5e9;
lamd  = 1.55e-6;
h     = 6.626e-34; cc = 3e8; q = 1.602e-19;
Rresp = 0.9;                                 % responsivity (A/W)
Id    = 5e-9;                                % dark current (A)
RL    = 50; T = 300; kB = 1.381e-23;
Fn    = 3;                                   % amplifier noise figure (linear)
Be    = Rb/2;

Prx_dBm = -40:0.05:-15;
Prx_W   = 1e-3*10.^(Prx_dBm/10);
Ip      = Rresp*Prx_W;
sig_sh  = sqrt(2*q*(Ip+Id)*Be);
sig_th  = sqrt(4*kB*T*Fn*Be/RL);
% Shot and thermal noise are statistically independent, so their variances
% add. The total rms noise is the square root of the sum of the squares,
% not the sum of the two rms values.
sig_tot = sqrt(sig_sh.^2 + sig_th.^2);
Qf      = Ip./sig_tot;
BER     = 0.5*erfc(Qf/sqrt(2));

fprintf('PART C : Digital link at %.1f Gbit/s\n',Rb/1e9);
fprintf('%14s %12s %16s\n','Power (dBm)','Q factor','BER');
for p = [-35 -30 -28 -26 -24 -22 -20]
    [~,i] = min(abs(Prx_dBm-p));
    fprintf('%14.1f %12.2f %16.3e\n',p,Qf(i),BER(i));
end
[~,i9] = min(abs(BER-1e-9));
fprintf('\nReceiver sensitivity at BER = 1e-9 : %.2f dBm (Q = %.2f)\n',...
        Prx_dBm(i9),Qf(i9));
fprintf('A bit error rate of 1e-9 corresponds to Q = 6, which is the\n');
fprintf('standard design point for a digital optical receiver.\n\n');

%% PART D : Analog link, carrier to noise ratio
m      = 0.5;                                % optical modulation index
RIN_dB = -155;                               % relative intensity noise (dB/Hz)
Ba     = 4e6;                                % analog bandwidth (Hz)
RIN    = 10^(RIN_dB/10);
Ipa    = Rresp*Prx_W;
Csig   = 0.5*(m*Ipa).^2;
Nrin   = RIN*Ipa.^2*Ba;
Nsh    = 2*q*Ipa*Ba;
Nth    = 4*kB*T*Fn*Ba/RL;
CNR    = 10*log10(Csig./(Nrin+Nsh+Nth));

fprintf('PART D : Analog link, modulation index %.1f\n',m);
fprintf('%14s %12s %12s %12s %12s\n',...
        'Power (dBm)','CNR (dB)','RIN share','Shot share','Thermal share');
for p = [-30 -25 -20 -15]
    [~,i] = min(abs(Prx_dBm-p));
    tot = Nrin(i)+Nsh(i)+Nth;
    fprintf('%14.1f %12.2f %12.3f %12.3f %12.3f\n',...
            p,CNR(i),Nrin(i)/tot,Nsh(i)/tot,Nth/tot);
end
fprintf('\nAt low received power the thermal noise of the receiver dominates\n');
fprintf('and the carrier to noise ratio rises 2 dB for every 1 dB of\n');
fprintf('optical power. At high power the relative intensity noise of the\n');
fprintf('source dominates and the carrier to noise ratio saturates.\n');

%% Figure 1 : attenuation spectrum
figure;
plot(lam,Atot,'LineWidth',2); hold on;
plot(lam,Aray,'--','LineWidth',1.4);
plot(lam,Air,'-.','LineWidth',1.4);
for k = 1:3, plot(win(k),alp(k),'ro','MarkerFaceColor','r','MarkerSize',6); end
grid on; axis([0.8 1.7 0 2]);
xlabel('Wavelength (\mum)'); ylabel('Attenuation (dB/km)');
legend('Total','Rayleigh scattering','Infrared absorption','Operating windows');
title('Optical Fibre Attenuation Spectrum');

%% Figure 2 : power budget against link length
figure;
Lk = 0:1:120;
for k = 1:3
    plot(Lk,Ptx - Nconn*Lconn - (alp(k)+Lsplice/Lspan)*Lk,'LineWidth',1.8); hold on;
end
plot([0 120],[Prx Prx],'k--','LineWidth',1.6);
plot([0 120],[Prx+Margin Prx+Margin],'k:','LineWidth',1.6);
grid on; axis([0 120 -45 5]);
xlabel('Link length (km)'); ylabel('Received power (dBm)');
legend('0.85 \mum','1.31 \mum','1.55 \mum','Sensitivity','Sensitivity + margin');
title('Optical Power Budget');

%% Figure 3 : digital link performance
figure;
semilogy(Prx_dBm,BER,'LineWidth',2); hold on;
semilogy([min(Prx_dBm) max(Prx_dBm)],[1e-9 1e-9],'k--','LineWidth',1.5);
plot([Prx_dBm(i9) Prx_dBm(i9)],[1e-14 1],'r:','LineWidth',1.5);
grid on; axis([-40 -15 1e-14 1]);
xlabel('Received optical power (dBm)'); ylabel('Bit error rate');
legend('BER','1 \times 10^{-9} reference','Receiver sensitivity');
title(sprintf('Digital Optical Link at %.1f Gbit/s',Rb/1e9));

%% Figure 4 : analog link carrier to noise ratio
figure;
plot(Prx_dBm,CNR,'LineWidth',2); hold on;
plot(Prx_dBm,10*log10(Csig./Nth),'--','LineWidth',1.4);
plot(Prx_dBm,10*log10(Csig./Nrin),'-.','LineWidth',1.4);
grid on; axis([-40 -15 0 60]);
xlabel('Received optical power (dBm)'); ylabel('CNR (dB)');
legend('Total CNR','Thermal noise limit','RIN limit','Location','southeast');
title('Analog Optical Link Carrier to Noise Ratio');
