% Experiment 2 : Free Space, Log-Distance and Two-Ray Ground Path Loss Models
clc; clear; close all;

%% Parameters
Pt   = 43;              % transmit power (dBm)
Gt   = 12;              % base station antenna gain (dBi)
Gr   = 0;               % mobile antenna gain (dBi)
fc   = 900e6;           % carrier frequency (Hz)
ht   = 30;              % base station antenna height (m)
hr   = 1.5;             % mobile antenna height (m)
d0   = 100;             % reference distance (m)
Sens = -104;            % receiver sensitivity (dBm)
c    = 3e8;
lambda = c/fc;

EIRP  = Pt + Gt + Gr;
maxPL = EIRP - Sens;

fprintf('Wavelength                    = %.4f m\n',lambda);
fprintf('Pt + Gt + Gr                  = %.1f dBm\n',EIRP);
fprintf('Maximum tolerable path loss = %.1f dB\n\n',maxPL);

d = logspace(1,5,4000);

%% Free space path loss
PL_free = 20*log10(4*pi*d/lambda);
PL_d0   = 20*log10(4*pi*d0/lambda);
fprintf('Free space path loss at 1 km = %.2f dB\n\n',20*log10(4*pi*1000/lambda));

%% Log-distance path loss
nvec   = [2 2.7 3.5 4];
PL_log = zeros(numel(nvec),numel(d));
for k = 1:numel(nvec)
    PL_log(k,:) = PL_d0 + 10*nvec(k)*log10(d/d0);
end

%% Two-ray ground reflection model
dc     = 4*pi*ht*hr/lambda;
dlos   = sqrt(d.^2 + (ht-hr)^2);
dref   = sqrt(d.^2 + (ht+hr)^2);
dphi   = 2*pi*(dref-dlos)/lambda;
Glin   = 10^((Gt+Gr)/10);
Efield = 1./dlos - (1./dref).*exp(-1j*dphi);
Pr_2ray = 10*log10(10^(Pt/10)*Glin*(lambda/(4*pi))^2 .* abs(Efield).^2);
Pr_asym = 10*log10(10^(Pt/10)*Glin*(ht^2*hr^2)./(d.^4));

fprintf('Two-ray crossover distance   = %.1f m (%.3f km)\n\n',dc,dc/1000);

%% Range at receiver sensitivity
fprintf('%-26s %10s %14s\n','Model','Exponent','Range (km)');
Rfree = d0*10^((maxPL - PL_d0)/20);
fprintf('%-26s %10.1f %14.2f\n','Free space',2.0,Rfree/1000);
for k = 1:numel(nvec)
    Rk = d0*10^((maxPL - PL_d0)/(10*nvec(k)));
    fprintf('%-26s %10.1f %14.2f\n',sprintf('Log-distance n = %.1f',nvec(k)),nvec(k),Rk/1000);
end
R2ray = 10^((EIRP + 20*log10(ht*hr) - Sens)/40);
fprintf('%-26s %10.1f %14.2f\n','Two-ray asymptotic',4.0,R2ray/1000);

%% Antenna height gain
fprintf('\nAntenna height gain 20log10(ht*hr) = %.2f dB\n',20*log10(ht*hr));
fprintf('Height 30 m to 45 m gains %.2f dB\n',20*log10(45/30));

%% Figure 1 : all models
Pr_free = EIRP - PL_free;
Pr_log  = EIRP - PL_log;
figure;
semilogx(d/1000,Pr_free,'LineWidth',2); hold on;
semilogx(d/1000,Pr_log(2,:),'LineWidth',1.6);
semilogx(d/1000,Pr_log(3,:),'LineWidth',1.6);
semilogx(d/1000,Pr_2ray,'LineWidth',1);
semilogx(d/1000,Pr_asym,'--','LineWidth',1.6);
semilogx([min(d) max(d)]/1000,[Sens Sens],'k:','LineWidth',1.8);
grid on; axis([0.01 100 -180 40]);
xlabel('Distance (km)'); ylabel('Received Power (dBm)');
legend('Free space','Log-distance n = 2.7','Log-distance n = 3.5',...
       'Two-ray exact','Two-ray asymptotic','Sensitivity','Location','southwest');
title('Received Power against Distance, 900 MHz');

%% Figure 2 : interference lobes below crossover
dz = linspace(50,4000,6000);
dl = sqrt(dz.^2+(ht-hr)^2); dr = sqrt(dz.^2+(ht+hr)^2);
Ez = 1./dl - (1./dr).*exp(-1j*2*pi*(dr-dl)/lambda);
Pz = 10*log10(10^(Pt/10)*Glin*(lambda/(4*pi))^2 .* abs(Ez).^2);
Pa = 10*log10(10^(Pt/10)*Glin*(ht^2*hr^2)./(dz.^4));
figure;
plot(dz/1000,Pz,'LineWidth',1.2); hold on;
plot(dz/1000,Pa,'--','LineWidth',1.8);
plot([dc dc]/1000,[min(Pz) max(Pz)],'k-.','LineWidth',1.4);
grid on; xlabel('Distance (km)'); ylabel('Received Power (dBm)');
legend('Two-ray exact','Asymptotic d^{-4}','Crossover distance');
title('Interference Lobes below the Crossover Distance');

%% Figure 3 : sensitivity of range to path loss exponent
nn = 2:0.05:5;
Rn = d0*10.^((maxPL - PL_d0)./(10*nn));
figure;
semilogy(nn,Rn/1000,'LineWidth',2); grid on;
xlabel('Path Loss Exponent n'); ylabel('Range (km)');
title('Range against Assumed Path Loss Exponent');
