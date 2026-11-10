% Experiment 9 : OFDM Transceiver - Orthogonality, Cyclic Prefix and PAPR
clc; clear; close all;
rand('state',9); randn('state',9);

%% Parameters
Nfft  = 64;             % FFT size
Ncp   = 16;             % cyclic prefix length
Nused = 52;             % used subcarriers
Nsym  = 2000;           % OFDM symbols
BW    = 20e6;           % signal bandwidth (Hz)
Mmod  = 4;              % QPSK

df  = BW/Nfft;
Tu  = 1/df;
Tcp = Ncp/Nfft*Tu;
Ts  = Tu + Tcp;

fprintf('FFT size                   = %d\n',Nfft);
fprintf('Used subcarriers           = %d\n',Nused);
fprintf('Subcarrier spacing         = %.4f kHz\n',df/1000);
fprintf('Useful symbol duration     = %.4f us\n',Tu*1e6);
fprintf('Cyclic prefix duration     = %.4f us\n',Tcp*1e6);
fprintf('Total symbol duration      = %.4f us\n',Ts*1e6);
fprintf('Cyclic prefix overhead     = %.2f %%\n',100*Ncp/(Nfft+Ncp));
fprintf('Symbol efficiency          = %.2f %%\n\n',100*Nfft/(Nfft+Ncp));

used = [Nfft/2-Nused/2+1 : Nfft/2, Nfft/2+2 : Nfft/2+1+Nused/2];

%% QPSK symbol generation
bits = randi([0 1],2,Nused*Nsym);
X    = ((1-2*bits(1,:)) + 1j*(1-2*bits(2,:)))/sqrt(2);
X    = reshape(X,Nused,Nsym);

%% OFDM modulation
grid_tx = zeros(Nfft,Nsym);
grid_tx(used,:) = X;
x  = ifft(ifftshift(grid_tx,1),Nfft)*sqrt(Nfft);
xc = [x(end-Ncp+1:end,:); x];
s  = xc(:);

%% Peak to average power ratio
papr = zeros(1,Nsym);
for k = 1:Nsym
    b = xc(:,k);
    papr(k) = max(abs(b).^2)/mean(abs(b).^2);
end
papr_dB = 10*log10(papr);
fprintf('PAPR mean                       = %.2f dB\n',mean(papr_dB));
fprintf('PAPR exceeded by 1 %% of symbols = %.2f dB\n',prctile_local(papr_dB,99));
fprintf('Theoretical limit 10log10(N) = %.2f dB\n\n',10*log10(Nfft));

%% Effect of cyclic prefix length on the error vector magnitude
% The metric is EVM rather than bit error rate. At the signal to noise
% ratio used here QPSK is robust enough that intersymbol interference from
% a short cyclic prefix degrades the constellation without yet producing
% bit errors, so EVM resolves the effect and BER does not.
h   = [1 0 0.6 0 0.3];                    % 5-tap multipath channel
Lch = numel(h);
cpv = [0 1 2 3 4 8 16];
snr = 30;
H   = fftshift(fft(h,Nfft)).';
evm = zeros(size(cpv));
fprintf('%10s %14s %14s %12s\n','CP length','CP >= L-1 ?','EVM','EVM (dB)');
for m = 1:numel(cpv)
    cp = cpv(m);
    if cp == 0, xcp = x; else, xcp = [x(end-cp+1:end,:); x]; end
    ss = xcp(:);
    rr = filter(h,1,ss);
    n  = sqrt(10^(-snr/10)*mean(abs(rr).^2)/2)*(randn(size(rr))+1j*randn(size(rr)));
    rr = rr + n;
    Lt = Nfft + cp;
    rm = reshape(rr(1:Nsym*Lt),Lt,Nsym);
    rm = rm(cp+1:end,:);
    Y  = fftshift(fft(rm,Nfft),1)/sqrt(Nfft);
    Ye = Y(used,:)./repmat(H(used),1,Nsym);
    evm(m) = sqrt(mean(abs(Ye(:)-X(:)).^2)/mean(abs(X(:)).^2));
    fprintf('%10d %14s %14.4f %12.2f\n',cp,...
            ternary_local(cp>=Lch-1,'yes','no'),evm(m),20*log10(evm(m)));
end
fprintf('\nThe channel delay spread is %d samples. Once the cyclic prefix\n',Lch-1);
fprintf('reaches that length the EVM stops improving, because orthogonality\n');
fprintf('has been restored and only thermal noise remains.\n\n');

%% Orthogonality check
fprintf('Subcarrier orthogonality check\n');
k1 = 10; k2 = 17;
v1 = exp(1j*2*pi*k1*(0:Nfft-1)/Nfft);
v2 = exp(1j*2*pi*k2*(0:Nfft-1)/Nfft);
fprintf('  Inner product of subcarriers %d and %d = %.3e\n',k1,k2,abs(v1*v2')/Nfft);
fprintf('  Inner product of subcarrier %d with itself = %.4f\n\n',k1,abs(v1*v1')/Nfft);

%% Figure 1 : PAPR complementary cumulative distribution
figure;
lv = 0:0.1:12;
cc = zeros(size(lv));
for k = 1:numel(lv), cc(k) = mean(papr_dB > lv(k)); end
semilogy(lv,max(cc,1e-4),'LineWidth',2); grid on;
axis([0 12 1e-4 1]);
xlabel('PAPR_0 (dB)'); ylabel('P ( PAPR > PAPR_0 )');
title(sprintf('PAPR Distribution, N = %d Subcarriers',Nfft));

%% Figure 2 : transmitted spectrum
nfft2 = 1024;
w = 0.5 - 0.5*cos(2*pi*(0:nfft2-1)'/(nfft2-1));
nseg = floor(numel(s)/nfft2); Ps = zeros(nfft2,1);
for m = 1:nseg
    Ps = Ps + abs(fft(s((m-1)*nfft2+1:m*nfft2).*w)).^2;
end
Ps = fftshift(Ps/nseg);
fa = (-nfft2/2:nfft2/2-1)'*BW/nfft2;
figure;
plot(fa/1e6,10*log10(Ps/max(Ps)),'LineWidth',1.2); grid on;
xlabel('Frequency (MHz)'); ylabel('Normalised PSD (dB)');
title('OFDM Transmitted Spectrum');

%% Figure 3 : effect of the cyclic prefix on the constellation
figure;
for idx = 1:2
    cpsel = [0 16]; cp = cpsel(idx);
    xcp = [x(end-max(cp,1)+1:end,:); x];
    if cp == 0, xcp = x; end
    ss = xcp(:); rr = filter(h,1,ss);
    n  = sqrt(10^(-snr/10)*mean(abs(rr).^2)/2)*(randn(size(rr))+1j*randn(size(rr)));
    rr = rr + n;
    Lt = Nfft + cp;
    rm = reshape(rr(1:Nsym*Lt),Lt,Nsym);
    rm = rm(cp+1:end,:);
    Y  = fftshift(fft(rm,Nfft),1)/sqrt(Nfft);
    Ye = Y(used,:)./repmat(H(used),1,Nsym);
    subplot(1,2,idx);
    plot(real(Ye(:)),imag(Ye(:)),'.','MarkerSize',1);
    axis square; axis([-2 2 -2 2]); grid on;
    xlabel('In-phase'); ylabel('Quadrature');
    title(sprintf('Cyclic prefix = %d samples',cp));
end

%% Figure 4 : EVM against cyclic prefix length
figure;
plot(cpv,20*log10(evm),'o-','LineWidth',2); hold on;
plot([Lch-1 Lch-1],[min(20*log10(evm))-2 max(20*log10(evm))+2],'k--','LineWidth',1.5);
grid on; xlabel('Cyclic prefix length (samples)'); ylabel('EVM (dB)');
legend('Measured EVM','Channel delay spread L-1','Location','northeast');
title('Cyclic Prefix Length against Error Vector Magnitude');

%% Figure 5 : time domain OFDM symbol with cyclic prefix
figure;
plot(abs(xc(:,1)),'LineWidth',1.4); hold on;
plot(1:Ncp,abs(xc(1:Ncp,1)),'r','LineWidth',2.2);
plot(Nfft+1:Nfft+Ncp,abs(xc(Nfft+1:Nfft+Ncp,1)),'r--','LineWidth',2.2);
grid on; xlabel('Sample'); ylabel('Magnitude');
legend('OFDM symbol','Cyclic prefix','Copied from the tail');
title('Cyclic Prefix Insertion');
