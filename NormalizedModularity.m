%NORMALIZED MODULARITY Generate Weighted Spectral Cluster Ensemble.
%
%   Input  : modules      : set modules as cell array of vectors, ex: {[1,2,3],[4,5,6]}
%            adj          : adjacency matrix
%
%   Output : normalized modularity metric, in [0,1]
%
% Written by: Muhammad Yousefnezhad,
%             iBRAIN group, Department of Computer Science & Technology,
%             Nanjing University of Aeronautics & Astronautics  
function NQ=NormalizedModularity(modules,adj)
nedges=numedges(adj); % total number of edges
Q = 0;
for m=1:length(modules)
  e_mm=numedges(adj(modules{m},modules{m}))/nedges;
  a_m=sum(sum(adj(modules{m},:)))/(2*nedges);
  Q = Q + (e_mm - a_m^2);
end
NQ = (Q + 1) / 2;
end


