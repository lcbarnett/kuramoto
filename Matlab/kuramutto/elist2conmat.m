function [C,n,k] = elist2conmat(elist)

% Convert edge list (cell array) to connectivity matrix

n = max(length(elist),max(cellfun(@max,elist))); % number of nodes
C = zeros(n);
for e = 1:n
	C(elist{e},e) = 1;
end
k = nnz(C); % number of edges
