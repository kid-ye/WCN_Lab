function Pb = erlangb(C,A)
% Blocking probability of a lost calls cleared system, stable recursion.
Pb = 1;
for k = 1:C
    Pb = (A*Pb)/(k + A*Pb);
end
end
