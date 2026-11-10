% Experiment 7 : GSM TDMA Frame Structure and Logical Channel Mapping
clc; clear; close all;

%% Parameters
Rb    = 270833.33;        % gross bit rate (bit/s)
Nslot = 8;                % timeslots per TDMA frame
Nbits = 148;              % bits in a normal burst
Guard = 8.25;             % guard period (bit periods)
Tb    = 1/Rb;

%% Burst composition
tail1 = 3; data1 = 57; steal1 = 1; tsc = 26; steal2 = 1; data2 = 57; tail2 = 3;
total = tail1+data1+steal1+tsc+steal2+data2+tail2;

fprintf('GSM normal burst composition\n');
fprintf('  Tail bits                 : %d + %d\n',tail1,tail2);
fprintf('  Encrypted data bits       : %d + %d\n',data1,data2);
fprintf('  Stealing flags            : %d + %d\n',steal1,steal2);
fprintf('  Training sequence         : %d\n',tsc);
fprintf('  Total                     : %d bits\n',total);
fprintf('  Guard period              : %.2f bit periods\n\n',Guard);

%% Timing hierarchy
Tbit     = Tb;
Tburst   = (Nbits + Guard)*Tb;
Tframe   = Nslot*Tburst;
Tmulti26 = 26*Tframe;
Tmulti51 = 51*Tframe;
Tsuper   = 51*Tmulti26;
Thyper   = 2048*Tsuper;

fprintf('GSM timing hierarchy\n');
fprintf('  Bit period                : %.4f us\n',Tbit*1e6);
fprintf('  Timeslot (burst) duration : %.4f us\n',Tburst*1e6);
fprintf('  TDMA frame (8 slots)      : %.4f ms\n',Tframe*1e3);
fprintf('  Traffic multiframe (26)   : %.2f ms\n',Tmulti26*1e3);
fprintf('  Control multiframe (51)   : %.2f ms\n',Tmulti51*1e3);
fprintf('  Superframe (51x26)        : %.3f s\n',Tsuper);
fprintf('  Hyperframe (2048 super)   : %.2f s\n\n',Thyper);

%% Useful bit rate of one traffic channel
useful = (data1+data2)*Nslot/Nslot;          % bits per burst carrying data
Rtch   = 114*24/(Tmulti26);                  % 24 traffic bursts per 26 frames
fprintf('Traffic channel gross rate  : %.2f kbit/s\n',114/Tframe/1000);
fprintf('Traffic channel useful rate : %.2f kbit/s\n\n',Rtch/1000);

%% 26-frame traffic multiframe map
tmap = repmat({'T'},1,26);
tmap{13} = 'S';           % SACCH
tmap{26} = 'I';           % idle
fprintf('26-frame traffic multiframe (T traffic, S SACCH, I idle)\n  ');
fprintf('%s ',tmap{:}); fprintf('\n\n');

%% 51-frame control multiframe map, downlink timeslot 0, non-combined
% Configuration: BCCH carrier, timeslot 0, downlink, non-combined, that is
% the SDCCH is carried on a separate timeslot. Frame numbers below are
% counted from 0 and may be checked directly against 3GPP TS 45.002.
cmap = cell(1,51);
for k = 1:51, cmap{k} = 'C'; end
cmap{1}  = 'F'; cmap{2}  = 'S';
cmap{11} = 'F'; cmap{12} = 'S';
cmap{21} = 'F'; cmap{22} = 'S';
cmap{31} = 'F'; cmap{32} = 'S';
cmap{41} = 'F'; cmap{42} = 'S';
cmap{51} = 'I';
for k = 3:6, cmap{k} = 'B'; end
fprintf('51-frame control multiframe, downlink timeslot 0, non-combined\n');
fprintf('  F FCCH   S SCH   B BCCH   C CCCH (PCH and AGCH)   I idle\n  ');
fprintf('%s ',cmap{:}); fprintf('\n');
lab = 'FSBCI';
for q = 1:numel(lab)
    idx = find(cellfun(@(z) z==lab(q), cmap)) - 1;
    fprintf('  %c at frames :',lab(q));
    fprintf(' %d',idx); fprintf('\n');
end
fprintf('  CCCH occupies %d of the 51 frames\n\n',sum(cellfun(@(z) z=='C',cmap)));

%% Logical channel summary
fprintf('%-8s %-12s %-10s %s\n','Channel','Direction','Type','Function');
fprintf('%-8s %-12s %-10s %s\n','FCCH','Downlink','Broadcast','Frequency correction');
fprintf('%-8s %-12s %-10s %s\n','SCH','Downlink','Broadcast','Frame synchronisation, BSIC');
fprintf('%-8s %-12s %-10s %s\n','BCCH','Downlink','Broadcast','Cell and network parameters');
fprintf('%-8s %-12s %-10s %s\n','PCH','Downlink','Common','Paging of mobile stations');
fprintf('%-8s %-12s %-10s %s\n','AGCH','Downlink','Common','Assignment of a dedicated channel');
fprintf('%-8s %-12s %-10s %s\n','RACH','Uplink','Common','Random access request');
fprintf('%-8s %-12s %-10s %s\n','SDCCH','Both','Dedicated','Call setup signalling');
fprintf('%-8s %-12s %-10s %s\n','SACCH','Both','Dedicated','Power control, timing advance');
fprintf('%-8s %-12s %-10s %s\n','TCH','Both','Traffic','Encoded speech or data');

%% Figure 1 : burst structure
figure;
seg = [tail1 data1 steal1 tsc steal2 data2 tail2 Guard];
lab = {'T','Data','S','Training','S','Data','T','Guard'};
col = [0.35 0.35 0.35;0.33 0.55 0.78;0.85 0.33 0.31;0.47 0.72 0.42;...
       0.85 0.33 0.31;0.33 0.55 0.78;0.35 0.35 0.35;0.85 0.85 0.85];
x = 0;
hold on;
for k = 1:numel(seg)
    rectangle('Position',[x 0 seg(k) 1],'FaceColor',col(k,:));
    if seg(k) > 3
        text(x+seg(k)/2,0.5,lab{k},'HorizontalAlignment','center','Color','w');
    end
    x = x + seg(k);
end
xlim([0 x]); ylim([0 1]); set(gca,'YTick',[]);
xlabel('Bit periods'); title('GSM Normal Burst, 156.25 Bit Periods');

%% Figure 2 : slot occupancy over one traffic multiframe
occ = zeros(Nslot,26);
for f = 1:26
    for s = 1:Nslot
        if strcmp(tmap{f},'T'), occ(s,f) = 1;
        elseif strcmp(tmap{f},'S'), occ(s,f) = 2;
        else, occ(s,f) = 0; end
    end
end
figure;
imagesc(1:26,1:Nslot,occ); set(gca,'YDir','normal');
colormap([0.9 0.9 0.9;0.33 0.55 0.78;0.95 0.75 0.30]);
xlabel('Frame number in the 26-frame multiframe'); ylabel('Timeslot');
title('Traffic Multiframe Occupancy (grey idle, blue TCH, orange SACCH)');

%% Figure 3 : timing hierarchy on a logarithmic scale
figure;
names = {'Bit','Timeslot','Frame','Multiframe 26','Multiframe 51','Superframe','Hyperframe'};
vals  = [Tbit Tburst Tframe Tmulti26 Tmulti51 Tsuper Thyper];
barh(log10(vals)); set(gca,'YTickLabel',names);
grid on; xlabel('log_{10} ( duration in seconds )');
title('GSM Timing Hierarchy');
