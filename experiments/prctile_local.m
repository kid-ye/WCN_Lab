function y = prctile_local(x,p)
x = sort(x(:));
n = numel(x);
y = interp1(100*((1:n)'-0.5)/n, x, min(max(p,100*0.5/n),100*(n-0.5)/n));
end
