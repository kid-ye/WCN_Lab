function Hm = bfs_hops(Adj)
% Hop count between every pair of nodes by breadth first search.
% Unreachable pairs are returned as Inf.
N  = size(Adj,1);
Hm = inf(N);
for s = 1:N
    d = inf(1,N); d(s) = 0; q = s;
    while ~isempty(q)
        v = q(1); q(1) = [];
        nb = find(Adj(v,:) & isinf(d));
        d(nb) = d(v) + 1;
        q = [q nb];
    end
    Hm(s,:) = d;
end
end
