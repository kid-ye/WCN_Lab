% Experiment 6 : Adaptive Equalisation of a Multipath Channel using the LMS Algorithm
clc; clear; close all;
randn('state',6);

%% Parameters
h      = [0.3 0.9 0.3];      % multipath channel impulse response
Meq    = 11;                 % equaliser length (taps)
D      = 6;                  % decision delay (samples)
mu     = 0.01;               % LMS step size
SNR_dB = 20;                 % signal to noise ratio
Nsym   = 5000;               % symbols per run
Nrun   = 200;                % runs averaged for the learning curve

Lh   = numel(h);
sigv2 = 10^(-SNR_dB/10)*sum(h.^2);

fprintf('Channel impulse response = [%s]\n',num2str(h));
fprintf('Equaliser length         = %d taps\n',Meq);
fprintf('Decision delay           = %d samples\n',D);
fprintf('Step size                = %.3f\n',mu);
fprintf('Noise variance           = %.5f\n\n',sigv2);

%% Wiener solution for comparison
H = zeros(Meq, Meq+Lh-1);
for k = 1:Meq
    H(k,k:k+Lh-1) = h;
end
R = H*H' + sigv2*eye(Meq);
p = H(:,D+1);
w_opt = R\p;
Jmin  = 1 - p'*w_opt;

fprintf('Minimum mean square error (Wiener) = %.6f\n',Jmin);
fprintf('Maximum eigenvalue of R            = %.4f\n',max(eig(R)));
fprintf('Step size stability limit          = %.4f\n\n',2/max(eig(R)));

%% LMS adaptation, averaged learning curve
J = zeros(1,Nsym);
for run = 1:Nrun
    s = 2*randi([0 1],1,Nsym+Meq+Lh) - 1;
    x = filter(h,1,s) + sqrt(sigv2)*randn(size(s));
    w = zeros(Meq,1);
    for n = Meq:Nsym
        u = x(n:-1:n-Meq+1).';
        y = w'*u;
        e = s(n-D) - y;
        w = w + mu*e*u;
        J(n) = J(n) + e^2;
    end
end
J = J/Nrun;

fprintf('Steady state mean square error (LMS) = %.6f\n',mean(J(end-499:end)));
fprintf('Tap error norm ||w - w_opt||         = %.6f\n\n',norm(w - w_opt));

%% Bit error rate against signal to noise ratio
% The channel peak occurs at tap index 1, so the unequalised detector is
% aligned to that tap. The equalised detector is aligned to the delay D.
[~,ipk]   = max(abs(h));  Dch = ipk - 1;
Ntest     = 2e5;
snr_sweep = 0:2:20;
ber_no    = zeros(size(snr_sweep));
ber_eq    = zeros(size(snr_sweep));
for k = 1:numel(snr_sweep)
    sv2 = 10^(-snr_sweep(k)/10)*sum(h.^2);
    s   = 2*randi([0 1],1,Ntest) - 1;
    x   = filter(h,1,s) + sqrt(sv2)*randn(1,Ntest);

    Rk  = H*H' + sv2*eye(Meq);
    wk  = Rk\p;
    y   = filter(wk,1,x);

    b   = (s+1)/2;
    ber_no(k) = mean(b(1:end-Dch) ~= (x(1+Dch:end) > 0));
    ber_eq(k) = mean(b(1:end-D)   ~= (y(1+D:end)   > 0));
end

fprintf('%10s %16s %16s\n','SNR (dB)','BER no equaliser','BER equalised');
for k = 1:numel(snr_sweep)
    fprintf('%10d %16.5f %16.5f\n',snr_sweep(k),ber_no(k),ber_eq(k));
end
fprintf('\n');

s = 2*randi([0 1],1,Ntest) - 1;
x = filter(h,1,s) + sqrt(sigv2)*randn(1,Ntest);
y = filter(w_opt,1,x);

%% Residual intersymbol interference
cimp = conv(h,w_opt);
[pk,ip] = max(abs(cimp));
isi_before = (sum(abs(h)) - max(abs(h)))/max(abs(h));
isi_after  = (sum(abs(cimp)) - pk)/pk;
fprintf('Intersymbol interference before equalisation = %.4f\n',isi_before);
fprintf('Intersymbol interference after equalisation  = %.4f\n',isi_after);

%% Figure 1 : learning curve
figure;
semilogy(J,'LineWidth',1.2); hold on;
semilogy([1 Nsym],[Jmin Jmin],'r--','LineWidth',1.8);
grid on; axis([1 Nsym 1e-2 2]);
xlabel('Iteration'); ylabel('Mean Square Error');
legend('LMS, averaged over 200 runs','Wiener minimum');
title('LMS Learning Curve');

%% Figure 2 : bit error rate against signal to noise ratio
figure;
semilogy(snr_sweep,ber_no,'o-','LineWidth',1.8); hold on;
semilogy(snr_sweep,max(ber_eq,1e-6),'s-','LineWidth',1.8);
semilogy(snr_sweep,0.5*erfc(sqrt(10.^(snr_sweep/10))),'k--','LineWidth',1.5);
grid on; axis([0 20 1e-6 1]);
xlabel('SNR (dB)'); ylabel('Bit Error Rate');
legend('Without equaliser','With equaliser','AWGN bound','Location','southwest');
title('Effect of Equalisation on Bit Error Rate');

%% Figure 3 : combined channel and equaliser response
figure;
subplot(1,2,1);
stem(0:Lh-1,h,'filled','LineWidth',1.5); grid on;
xlabel('Tap'); ylabel('Amplitude'); title('Channel Impulse Response');
subplot(1,2,2);
stem(0:numel(cimp)-1,cimp,'filled','LineWidth',1.5); grid on;
xlabel('Tap'); ylabel('Amplitude'); title('Channel and Equaliser Combined');

%% Figure 4 : signal at three points
figure;
subplot(3,1,1); plot(s(1:200),'LineWidth',1.2); grid on;
ylabel('Amplitude'); title('Transmitted BPSK Symbols');
subplot(3,1,2); plot(x(1:200),'LineWidth',1.2); grid on;
ylabel('Amplitude'); title('After the Multipath Channel');
subplot(3,1,3); plot(y(1+D:200+D),'LineWidth',1.2); grid on;
xlabel('Symbol Index'); ylabel('Amplitude'); title('Equaliser Output');

%% Figure 5 : equaliser tap weights
figure;
stem(0:Meq-1,w_opt,'filled','LineWidth',1.5); hold on;
stem(0:Meq-1,w,'r','LineWidth',1);
grid on; xlabel('Tap'); ylabel('Weight');
legend('Wiener solution','LMS after adaptation');
title('Equaliser Tap Weights');
