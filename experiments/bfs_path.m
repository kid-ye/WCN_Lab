function [p,h] = bfs_path(Adj,s,t)
% Shortest hop path from s to t by breadth first search.
N = size(Adj,1);
prev = zeros(1,N); d = inf(1,N); d(s) = 0; q = s;
while ~isempty(q)
    v = q(1); q(1) = [];
    nb = find(Adj(v,:) & isinf(d));
    d(nb) = d(v) + 1; prev(nb) = v;
    q = [q nb];
end
if isinf(d(t)), p = []; h = Inf; return; end
p = t; while p(1) ~= s, p = [prev(p(1)) p]; end
h = d(t);
end
