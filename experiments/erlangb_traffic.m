function A = erlangb_traffic(C,Pt)
% Offered traffic carried by C channels at blocking probability Pt.
if C == 0, A = 0; return; end
lo = 0; hi = 2*C;
while erlangb(C,hi) < Pt, hi = 2*hi; end
for k = 1:200
    mid = 0.5*(lo+hi);
    if erlangb(C,mid) < Pt, lo = mid; else hi = mid; end
end
A = 0.5*(lo+hi);
end
