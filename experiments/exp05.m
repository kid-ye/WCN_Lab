% Experiment 5 : Diversity Techniques - Selection, Equal Gain and Maximal Ratio Combining
clc; clear; close all;
randn('state',5);

%% Parameters
Lmax    = 4;                     % maximum number of diversity branches
Nsym    = 1e6;                   % channel realisations
gbar_dB = 0:2:30;                % average branch SNR (dB)
gth_dB  = 5;                     % outage threshold (dB)
gth     = 10^(gth_dB/10);

fprintf('Outage threshold = %d dB\n',gth_dB);
fprintf('Branch fading    = Rayleigh, branches independent\n\n');

Pout_sc  = zeros(Lmax,numel(gbar_dB));
Pout_egc = Pout_sc;  Pout_mrc = Pout_sc;
Tout_sc  = Pout_sc;  Tout_mrc = Pout_sc;

for iL = 1:Lmax
    L = iL;
    for ig = 1:numel(gbar_dB)
        gbar = 10^(gbar_dB(ig)/10);

        % Rayleigh branch amplitudes, unit mean power
        hmat = (randn(L,Nsym) + 1j*randn(L,Nsym))/sqrt(2);
        g    = gbar*abs(hmat).^2;                       % branch SNRs

        % Selection combining : strongest branch
        g_sc  = max(g,[],1);

        % Maximal ratio combining : SNRs add
        g_mrc = sum(g,1);

        % Equal gain combining : amplitudes add coherently
        g_egc = (sum(abs(hmat),1).^2)*gbar/L;

        Pout_sc(iL,ig)  = mean(g_sc  < gth);
        Pout_egc(iL,ig) = mean(g_egc < gth);
        Pout_mrc(iL,ig) = mean(g_mrc < gth);

        % Closed form results
        x = gth/gbar;
        Tout_sc(iL,ig)  = (1 - exp(-x))^L;
        s = 0;
        for k = 0:L-1, s = s + x^k/factorial(k); end
        Tout_mrc(iL,ig) = 1 - exp(-x)*s;
    end
end

%% Monte Carlo resolution
fprintf('Realisations per point = %g, resolution = %.1e\n\n',Nsym,1/Nsym);

%% Outage probability at a representative average SNR
idx = find(gbar_dB == 14);
fprintf('Outage probability at average branch SNR of %d dB\n',gbar_dB(idx));
fprintf('%10s %12s %12s %12s %12s %12s\n','Branches','SC sim','SC theory','EGC sim','MRC sim','MRC theory');
for L = 1:Lmax
    fprintf('%10d %12.3e %12.3e %12.3e %12.3e %12.3e\n',L,...
        Pout_sc(L,idx),Tout_sc(L,idx),Pout_egc(L,idx),Pout_mrc(L,idx),Tout_mrc(L,idx));
end

%% Diversity order from the slope of the outage curve
% Taken from the closed form results, which remain accurate below the
% Monte Carlo resolution of 1/Nsym.
fprintf('\nDiversity order, magnitude of the outage slope per decade of SNR\n');
fprintf('%10s %12s %12s %12s\n','Branches','SC','EGC','MRC');
i1 = find(gbar_dB == 20); i2 = find(gbar_dB == 30);
j1 = find(gbar_dB == 8);  j2 = find(gbar_dB == 14);
for L = 1:Lmax
    slope_sc  = log10(Tout_sc(L,i1))  - log10(Tout_sc(L,i2));
    slope_mrc = log10(Tout_mrc(L,i1)) - log10(Tout_mrc(L,i2));
    % EGC has no closed form, so its slope is taken from simulated points
    % that lie well above the Monte Carlo resolution.
    slope_egc = (log10(Pout_egc(L,j1)) - log10(Pout_egc(L,j2))) / ...
                ((gbar_dB(j2)-gbar_dB(j1))/10);
    fprintf('%10d %12.2f %12.2f %12.2f\n',L,slope_sc,slope_egc,slope_mrc);
end

%% Array gain : average combined SNR relative to one branch
fprintf('\nArray gain, average output SNR relative to a single branch (dB)\n');
fprintf('%9s %9s %9s %9s %9s %9s %9s\n','Branches',...
        'SC sim','SC thy','EGC sim','EGC thy','MRC sim','MRC thy');
gbar = 1;
for L = 1:Lmax
    hmat  = (randn(L,Nsym) + 1j*randn(L,Nsym))/sqrt(2);
    g     = gbar*abs(hmat).^2;
    a_sc  = 10*log10(mean(max(g,[],1))/gbar);
    a_egc = 10*log10(mean(sum(abs(hmat),1).^2)/(L*gbar));
    a_mrc = 10*log10(mean(sum(g,1))/gbar);
    t_sc  = 10*log10(sum(1./(1:L)));
    t_egc = 10*log10(1 + (L-1)*pi/4);
    t_mrc = 10*log10(L);
    fprintf('%9d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n',...
            L,a_sc,t_sc,a_egc,t_egc,a_mrc,t_mrc);
end

% Note : at very low outage probabilities a finite Monte Carlo simulation
% may record no outage events. Such points are not zero probability, they
% are simply below the resolution of 1/Nsym, and are excluded from the plot.

%% Replace unresolved zeros so that the logarithmic plots are valid
Pout_sc(Pout_sc==0)   = NaN;
Pout_egc(Pout_egc==0) = NaN;
Pout_mrc(Pout_mrc==0) = NaN;

%% Figure 1 : outage probability
figure;
mk = {'o-','s-','d-','^-'};
for L = 1:Lmax
    semilogy(gbar_dB,Pout_mrc(L,:),mk{L},'LineWidth',1.6); hold on;
end
semilogy(gbar_dB,Tout_mrc(1,:),'k--','LineWidth',1.4);
for L = 2:Lmax
    semilogy(gbar_dB,Tout_mrc(L,:),'k--','LineWidth',1.4);
end
grid on; axis([0 30 1e-6 1]);
xlabel('Average Branch SNR (dB)'); ylabel('Outage Probability');
legend('MRC L = 1','MRC L = 2','MRC L = 3','MRC L = 4','Theory','Location','southwest');
title('Maximal Ratio Combining : Outage Probability');

%% Figure 2 : comparison of the three combiners at L = 3
figure;
semilogy(gbar_dB,Pout_sc(1,:),'k-','LineWidth',2); hold on;
semilogy(gbar_dB,Pout_sc(3,:),'o-','LineWidth',1.6);
semilogy(gbar_dB,Pout_egc(3,:),'s-','LineWidth',1.6);
semilogy(gbar_dB,Pout_mrc(3,:),'d-','LineWidth',1.6);
grid on; axis([0 30 1e-6 1]);
xlabel('Average Branch SNR (dB)'); ylabel('Outage Probability');
legend('No diversity','Selection L = 3','Equal gain L = 3','Maximal ratio L = 3',...
       'Location','southwest');
title('Comparison of Combining Techniques');
