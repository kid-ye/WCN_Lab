% Experiment 3 : Okumura-Hata Model and Coverage Prediction
clc; clear; close all;

%% Parameters
f     = 900;            % carrier frequency (MHz)
hb    = 30;             % base station antenna height (m)
hm    = 1.5;            % mobile antenna height (m)
EIRP  = 55;             % Pt + Gt + Gr (dBm)
Sens  = -104;           % receiver sensitivity (dBm)
side  = 20;             % service area side (km)
Area  = side^2;
maxPL = EIRP - Sens;

d = 1:0.05:20;

%% Mobile antenna correction factor, medium city
a_hm = (1.1*log10(f) - 0.7)*hm - (1.56*log10(f) - 0.8);

%% Hata path loss
L_urban = 69.55 + 26.16*log10(f) - 13.82*log10(hb) - a_hm ...
          + (44.9 - 6.55*log10(hb))*log10(d);
L_sub   = L_urban - 2*(log10(f/28))^2 - 5.4;
L_open  = L_urban - 4.78*(log10(f))^2 + 18.33*log10(f) - 40.94;

nEff = (44.9 - 6.55*log10(hb))/10;

fprintf('Mobile antenna correction a(hm) = %.4f dB\n',a_hm);
fprintf('Effective path loss exponent    = %.4f\n',nEff);
fprintf('Maximum tolerable path loss     = %.1f dB\n\n',maxPL);

fprintf('%-12s %12s %12s %14s %14s\n','Environment','L(1 km)','L(10 km)','Radius (km)','Sites');
L1 = [L_urban(1) L_sub(1) L_open(1)];
env = {'Urban','Suburban','Open'};
Rcov = zeros(1,3); Nsite = zeros(1,3);
for k = 1:3
    Rcov(k)  = 10^((maxPL - L1(k))/(10*nEff));
    Acell    = 1.5*sqrt(3)*Rcov(k)^2;
    Nsite(k) = max(1,ceil(Area/Acell));
    L10 = L1(k) + 10*nEff;
    fprintf('%-12s %12.2f %12.2f %14.2f %14d\n',env{k},L1(k),L10,Rcov(k),Nsite(k));
end

fprintf('\nSuburban correction = %.2f dB\n',L_urban(1)-L_sub(1));
fprintf('Open area correction = %.2f dB\n',L_urban(1)-L_open(1));

fprintf('\nValidity of the Hata model : 150-1500 MHz, hb 30-200 m, hm 1-10 m, d 1-20 km\n');
for k = 1:3
    if Rcov(k) > 20
        fprintf('  %-9s radius %.1f km exceeds the 20 km validity limit\n',env{k},Rcov(k));
    end
end

%% Effect of base station antenna height
fprintf('\n%8s %10s %14s %10s\n','hb (m)','n','Radius (km)','Sites');
for h = [30 45 60 100]
    ne  = (44.9 - 6.55*log10(h))/10;
    Lu1 = 69.55 + 26.16*log10(f) - 13.82*log10(h) ...
          - ((1.1*log10(f)-0.7)*hm - (1.56*log10(f)-0.8));
    Ls1 = Lu1 - 2*(log10(f/28))^2 - 5.4;
    Rh  = 10^((maxPL - Ls1)/(10*ne));
    fprintf('%8d %10.4f %14.2f %10d\n',h,ne,Rh,max(1,ceil(Area/(1.5*sqrt(3)*Rh^2))));
end

%% Figure 1 : path loss by environment
figure;
plot(d,L_urban,'LineWidth',2); hold on;
plot(d,L_sub,'LineWidth',2);
plot(d,L_open,'LineWidth',2);
plot([1 20],[maxPL maxPL],'k--','LineWidth',1.8);
grid on; xlabel('Distance (km)'); ylabel('Median Path Loss (dB)');
legend('Urban','Suburban','Open area','Maximum tolerable','Location','southeast');
title('Okumura-Hata Median Path Loss, 900 MHz, h_b = 30 m');

%% Figure 2 : coverage map built from the calculated coverage radius
% The layout is generated for the URBAN environment, whose coverage radius
% is the smallest of the three and therefore the demanding case. Base
% station spacing is sqrt(3)*R so that the hexagons of circumradius R
% tessellate without gaps.
Rcell = Rcov(1);                          % urban coverage radius (km)
q = []; r = [];
nring = ceil(side/(sqrt(3)*Rcell)) + 1;
for a = -nring:nring
    for b = -nring:nring
        if abs(a+b) <= nring, q(end+1,1) = a; r(end+1,1) = b; end
    end
end
xb = sqrt(3)*Rcell*(q + r/2);  yb = 1.5*Rcell*r;
keep = abs(xb) <= side/2 + Rcell & abs(yb) <= side/2 + Rcell;
xb = xb(keep); yb = yb(keep);

g  = -side/2:0.1:side/2;
[GX,GY] = meshgrid(g,g);
Prx = -inf(size(GX));
for b = 1:numel(xb)
    dd  = max(sqrt((GX-xb(b)).^2 + (GY-yb(b)).^2), 1);
    Lb  = L_urban(1) + 10*nEff*log10(dd);
    Prx = max(Prx, EIRP - Lb);
end
covered = 100*mean(Prx(:) >= Sens);

figure;
imagesc(g,g,Prx); set(gca,'YDir','normal'); hold on;
plot(xb,yb,'kv','MarkerFaceColor','k','MarkerSize',7);
colorbar; axis equal; axis([-side/2 side/2 -side/2 side/2]);
xlabel('Distance (km)'); ylabel('Distance (km)');
title(sprintf('Urban Coverage : %d sites, R = %.2f km, %.1f %% above sensitivity',...
      numel(xb),Rcell,covered));

fprintf('\nCoverage map, urban environment\n');
fprintf('  Coverage radius used      = %.2f km\n',Rcell);
fprintf('  Base station spacing      = %.2f km\n',sqrt(3)*Rcell);
fprintf('  Sites placed in the area  = %d\n',numel(xb));
fprintf('  Area above sensitivity    = %.2f %%\n',covered);

%% Estimated against actual site count
Aest = ceil(Area/(1.5*sqrt(3)*Rcell^2));
fprintf('\n  Site count from the area ratio = %d\n',Aest);
fprintf('  Sites actually placed          = %d\n',numel(xb));
fprintf('  The area ratio underestimates the requirement because the hexagons\n');
fprintf('  at the boundary of a square service area are only partly used.\n');
