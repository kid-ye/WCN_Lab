% Experiment 11 : LTE Resource Grid, Pilot Insertion and Channel Estimation
clc; clear; close all;
randn('state',11); rand('state',11);

%% Parameters
Nsc    = 12;           % subcarriers per resource block
Nsym   = 7;            % OFDM symbols per slot, normal cyclic prefix
Nrb    = 25;           % resource blocks (5 MHz LTE)
df     = 15e3;         % subcarrier spacing (Hz)
Nfft   = 512;
pilotK = 6;            % pilot spacing in frequency (subcarriers)
pilotL = [1 5];        % pilot positions in time (symbol index)
snr_dB = 0:3:30;
Nslot  = 200;

Ntot  = Nrb*Nsc;
BW    = Ntot*df;
Tslot = 0.5e-3;

fprintf('Subcarrier spacing          = %.1f kHz\n',df/1000);
fprintf('Resource blocks             = %d\n',Nrb);
fprintf('Occupied subcarriers        = %d\n',Ntot);
fprintf('Occupied bandwidth          = %.3f MHz\n',BW/1e6);
fprintf('Symbols per slot            = %d\n',Nsym);
fprintf('Slot duration               = %.2f ms\n',Tslot*1e3);
fprintf('Resource elements per slot = %d\n',Ntot*Nsym);
fprintf('Pilot spacing in frequency = %d subcarriers\n',pilotK);
fprintf('Pilot symbols per slot      = %d of %d\n',numel(pilotL),Nsym);
pilotFrac = (numel(pilotL)*ceil(Ntot/pilotK))/(Ntot*Nsym);
fprintf('Pilot overhead              = %.2f %%\n\n',100*pilotFrac);

%% Coherence bandwidth check
tauRMS = 0.5e-6;                         % rms delay spread (s)
Bc     = 1/(5*tauRMS);
fprintf('Assumed rms delay spread   = %.2f us\n',tauRMS*1e6);
fprintf('Coherence bandwidth        = %.1f kHz\n',Bc/1000);
fprintf('Pilot spacing in frequency = %.1f kHz\n',pilotK*df/1000);
if pilotK*df < Bc
    fprintf('Pilot spacing is inside the coherence bandwidth, so linear\n');
    fprintf('interpolation between pilots is valid.\n\n');
else
    fprintf('Pilot spacing exceeds the coherence bandwidth. Interpolation\n');
    fprintf('will fail and the estimate will be biased.\n\n');
end

%% Channel estimation performance
pidx = 1:pilotK:Ntot;
mse  = zeros(size(snr_dB));
evm  = zeros(size(snr_dB));
for k = 1:numel(snr_dB)
    snr = 10^(snr_dB(k)/10);
    e   = 0; ev = 0;
    for m = 1:Nslot
        hd = (randn(1,6)+1j*randn(1,6)).*exp(-(0:5)/3)/sqrt(2);
        H  = fft(hd,Nfft); H = H(1:Ntot);
        Hg = repmat(H.',1,Nsym);

        b  = randi([0 1],2,Ntot*Nsym);
        X  = reshape(((1-2*b(1,:))+1j*(1-2*b(2,:)))/sqrt(2),Ntot,Nsym);
        X(pidx,pilotL) = (1+1j)/sqrt(2);

        N  = (randn(Ntot,Nsym)+1j*randn(Ntot,Nsym))/sqrt(2*snr);
        Y  = Hg.*X + N;

        Hp   = Y(pidx,pilotL)./X(pidx,pilotL);
        Hp   = mean(Hp,2);
        He   = interp1(pidx,Hp,1:Ntot,'linear','extrap').';
        Hest = repmat(He,1,Nsym);

        e  = e + mean(abs(Hest(:)-Hg(:)).^2)/mean(abs(Hg(:)).^2);
        % MMSE equalisation. Zero forcing divides by the channel and so
        % amplifies noise without limit on subcarriers in a deep fade.
        Xe = Y.*conj(Hest)./(abs(Hest).^2 + 1/snr);
        dm = setdiff(1:Nsym,pilotL);
        ev = ev + mean(abs(Xe(:,dm)-X(:,dm)).^2,'all')/mean(abs(X(:,dm)).^2,'all');
    end
    mse(k) = e/Nslot; evm(k) = sqrt(ev/Nslot);
end

fprintf('%10s %16s %14s %12s\n','SNR (dB)','Estimation MSE','EVM','EVM (dB)');
for k = 1:numel(snr_dB)
    fprintf('%10d %16.5f %14.4f %12.2f\n',snr_dB(k),mse(k),evm(k),20*log10(evm(k)));
end
fprintf('\nThe estimation error falls with signal to noise ratio until it\n');
fprintf('reaches the floor set by interpolation between pilots.\n');
fprintf('MMSE equalisation is used, since zero forcing amplifies noise\n');
fprintf('without limit on subcarriers lying in a deep fade.\n\n');

%% Figure 1 : resource grid with pilots
figure;
grid_map = zeros(Nsc*2,Nsym);
grid_map(:) = 1;
for c = pilotL
    grid_map(1:pilotK:end,c) = 2;
end
imagesc(1:Nsym,1:Nsc*2,grid_map); set(gca,'YDir','normal');
colormap([0.80 0.86 0.92;0.90 0.55 0.20]);
set(gca,'XTick',1:Nsym,'YTick',1:2:Nsc*2);
xlabel('OFDM symbol index'); ylabel('Subcarrier index');
title('LTE Resource Grid : Two Resource Blocks with Reference Signals');
grid on;

%% Figure 2 : channel estimation error
figure;
semilogy(snr_dB,mse,'o-','LineWidth',1.8); hold on;
semilogy(snr_dB,10.^(-snr_dB/10),'k--','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Normalised Mean Square Error');
legend('Least squares estimate with interpolation','1 / SNR reference');
title('Channel Estimation Error');

%% Figure 3 : error vector magnitude
figure;
plot(snr_dB,20*log10(evm),'s-','LineWidth',1.8); grid on;
xlabel('SNR (dB)'); ylabel('EVM (dB)');
title('Error Vector Magnitude after Equalisation');

%% Figure 4 : one channel realisation against its estimate
hd = (randn(1,6)+1j*randn(1,6)).*exp(-(0:5)/3)/sqrt(2);
H  = fft(hd,Nfft); H = H(1:Ntot).';
snr = 10^(20/10);
Xp = (1+1j)/sqrt(2)*ones(numel(pidx),1);
Yp = H(pidx).*Xp + (randn(numel(pidx),1)+1j*randn(numel(pidx),1))/sqrt(2*snr);
Hp = Yp./Xp;
He = interp1(pidx,Hp,1:Ntot,'linear','extrap').';
figure;
plot(20*log10(abs(H)),'LineWidth',1.6); hold on;
plot(20*log10(abs(He)),'--','LineWidth',1.6);
plot(pidx,20*log10(abs(Hp)),'o','MarkerSize',5);
grid on; xlabel('Subcarrier index'); ylabel('|H| (dB)');
legend('True channel','Interpolated estimate','Pilot measurements');
title('Frequency Selective Channel and its Estimate at 20 dB');
