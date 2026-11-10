% Experiment 10 : MIMO - SISO, Alamouti Transmit Diversity and Spatial Multiplexing
clc; clear; close all;
randn('state',10);

%% Parameters
snr_dB = 0:2:24;
Nsym   = 1e6;
Ncap   = 2e4;

fprintf('Rayleigh flat fading, perfect channel knowledge at the receiver\n');
fprintf('Modulation : BPSK for the error rate study\n');
fprintf('Symbols per point = %g, Monte Carlo resolution = %.1e\n\n',Nsym,1/Nsym);

ber_siso = zeros(size(snr_dB));
ber_ala  = zeros(size(snr_dB));
ber_mrc  = zeros(size(snr_dB));

for k = 1:numel(snr_dB)
    snr = 10^(snr_dB(k)/10);
    b   = randi([0 1],1,Nsym);
    s   = 1 - 2*b;

    % 1 x 1 SISO
    h = (randn(1,Nsym) + 1j*randn(1,Nsym))/sqrt(2);
    n = (randn(1,Nsym) + 1j*randn(1,Nsym))/sqrt(2*snr);
    y = h.*s + n;
    ber_siso(k) = mean(b ~= (real(conj(h).*y) < 0));

    % 2 x 1 Alamouti, transmit power split between the two antennas
    Np = Nsym/2;
    s1 = s(1:2:end); s2 = s(2:2:end);
    h1 = (randn(1,Np) + 1j*randn(1,Np))/sqrt(2);
    h2 = (randn(1,Np) + 1j*randn(1,Np))/sqrt(2);
    n1 = (randn(1,Np) + 1j*randn(1,Np))/sqrt(2*snr);
    n2 = (randn(1,Np) + 1j*randn(1,Np))/sqrt(2*snr);
    y1 = (h1.*s1 + h2.*s2)/sqrt(2) + n1;
    y2 = (-h1.*conj(s2) + h2.*conj(s1))/sqrt(2) + n2;
    z1 = conj(h1).*y1 + h2.*conj(y2);
    z2 = conj(h2).*y1 - h1.*conj(y2);
    bh = zeros(1,Nsym);
    bh(1:2:end) = real(z1) < 0;
    bh(2:2:end) = real(z2) < 0;
    ber_ala(k) = mean(b ~= bh);

    % 1 x 2 receive diversity with maximal ratio combining
    ha = (randn(1,Nsym) + 1j*randn(1,Nsym))/sqrt(2);
    hb = (randn(1,Nsym) + 1j*randn(1,Nsym))/sqrt(2);
    na = (randn(1,Nsym) + 1j*randn(1,Nsym))/sqrt(2*snr);
    nb = (randn(1,Nsym) + 1j*randn(1,Nsym))/sqrt(2*snr);
    z  = conj(ha).*(ha.*s + na) + conj(hb).*(hb.*s + nb);
    ber_mrc(k) = mean(b ~= (real(z) < 0));
end

fprintf('%10s %14s %14s %14s\n','SNR (dB)','SISO','Alamouti 2x1','MRC 1x2');
for k = 1:numel(snr_dB)
    fprintf('%10d %14.5f %14.5f %14.5f\n',snr_dB(k),ber_siso(k),ber_ala(k),ber_mrc(k));
end

%% Diversity order from the slope at high signal to noise ratio
% The slope is measured between 10 and 16 dB, where all three curves lie
% well above the Monte Carlo resolution of 1/Nsym.
i1 = find(snr_dB == 10); i2 = find(snr_dB == 16);
d_siso = (log10(ber_siso(i1))-log10(ber_siso(i2)))/((snr_dB(i2)-snr_dB(i1))/10);
d_ala  = (log10(ber_ala(i1)) -log10(ber_ala(i2)) )/((snr_dB(i2)-snr_dB(i1))/10);
d_mrc  = (log10(ber_mrc(i1)) -log10(ber_mrc(i2)) )/((snr_dB(i2)-snr_dB(i1))/10);
fprintf('\nDiversity order from the slope between 10 and 16 dB\n');
fprintf('  SISO          : %.2f  (expected 1)\n',d_siso);
fprintf('  Alamouti 2x1  : %.2f  (expected 2)\n',d_ala);
fprintf('  MRC 1x2       : %.2f  (expected 2)\n',d_mrc);
fprintf('\nAlamouti achieves the same diversity order as two-branch receive\n');
fprintf('combining but is 3 dB worse, because the transmit power is divided\n');
fprintf('between the two transmit antennas.\n\n');

%% Ergodic capacity
fprintf('%10s %12s %12s %12s %12s\n','SNR (dB)','1x1','1x2','2x1','2x2');
cap = zeros(4,numel(snr_dB));
for k = 1:numel(snr_dB)
    snr = 10^(snr_dB(k)/10);
    c11 = 0; c12 = 0; c21 = 0; c22 = 0;
    for m = 1:Ncap
        H11 = (randn + 1j*randn)/sqrt(2);
        c11 = c11 + log2(1 + snr*abs(H11)^2);
        H12 = (randn(2,1) + 1j*randn(2,1))/sqrt(2);
        c12 = c12 + log2(1 + snr*(H12'*H12));
        H21 = (randn(1,2) + 1j*randn(1,2))/sqrt(2);
        c21 = c21 + log2(1 + snr/2*(H21*H21'));
        H22 = (randn(2,2) + 1j*randn(2,2))/sqrt(2);
        c22 = c22 + real(log2(det(eye(2) + snr/2*(H22*H22'))));
    end
    cap(:,k) = [c11 c12 c21 c22]'/Ncap;
    fprintf('%10d %12.3f %12.3f %12.3f %12.3f\n',snr_dB(k),cap(1,k),cap(2,k),cap(3,k),cap(4,k));
end
fprintf('\nAt high signal to noise ratio the 2x2 capacity grows at twice the\n');
fprintf('rate of the 1x1 capacity, which is the spatial multiplexing gain.\n');
sl = (cap(:,end)-cap(:,end-1))/((snr_dB(end)-snr_dB(end-1))/3.0103);
fprintf('%10s %12s %12s %12s %12s\n','Slope','1x1','1x2','2x1','2x2');
fprintf('%10s %12.2f %12.2f %12.2f %12.2f\n','bit/s/Hz per 3 dB',sl(1),sl(2),sl(3),sl(4));

%% Figure 1 : bit error rate
figure;
semilogy(snr_dB,ber_siso,'o-','LineWidth',1.8); hold on;
semilogy(snr_dB,ber_ala,'s-','LineWidth',1.8);
semilogy(snr_dB,ber_mrc,'d-','LineWidth',1.8);
gr = 10.^(snr_dB/10);
semilogy(snr_dB,0.5*(1-sqrt(gr./(1+gr))),'k:','LineWidth',1.6);
semilogy(snr_dB,0.5*erfc(sqrt(gr)),'k--','LineWidth',1.4);
grid on; axis([0 24 1e-6 1]);
xlabel('SNR (dB)'); ylabel('Bit Error Rate');
legend('SISO 1x1','Alamouti 2x1','MRC 1x2','SISO theory','AWGN, no fading',...
       'Location','southwest');
title('Effect of Transmit and Receive Diversity');

%% Figure 2 : ergodic capacity
figure;
plot(snr_dB,cap(1,:),'o-','LineWidth',1.8); hold on;
plot(snr_dB,cap(2,:),'s-','LineWidth',1.8);
plot(snr_dB,cap(3,:),'d-','LineWidth',1.8);
plot(snr_dB,cap(4,:),'^-','LineWidth',1.8);
grid on; xlabel('SNR (dB)'); ylabel('Ergodic Capacity (bit/s/Hz)');
legend('1 x 1','1 x 2','2 x 1','2 x 2','Location','northwest');
title('Ergodic Capacity of Rayleigh MIMO Channels');
