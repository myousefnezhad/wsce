%WSCE Generate Weighted Spectral Cluster Ensemble.
%
%   Input  : data         : N-by-D data matrix, where N is the number of data,
%                           D is the number of dimensions
%            k            : number of clusters in the final result 
%            num_neighbors: number of nearest neighbors
%                           by default is 10     
%            block_size   : block size for partitioning data matrix
%                           by default is 10     
%            showden      : 1 for show dendrogram of the final result
%                           by default is 0     
%
%   Output : final result of clustering
% Written by: Muhammad Yousefnezhad,
%             iBRAIN group, Department of Computer Science & Technology,
%             Nanjing University of Aeronautics & Astronautics  
function [FinalResult] = WSCE(data, k, block_size, num_neighbors, showden) 
%% Initialize inputs
switch nargin
    case 2
        block_size = 10;
        num_neighbors = 10;
        showden = 0;
    case 3
        num_neighbors = 10;
        showden = 0;
    case 4
        showden = 0;
end
%% Generate (t-nearest-neighbor) sparse distance matrix.
Adj = gen_nn_distance(data,num_neighbors,block_size);
Adj = full(Adj);
%% Generate individul clusering results
WC = 0;
for (i = 2:k)
    [Ind, Mod] = wisedSpectral(Adj, i, 30);
    NQ = NormalizedModularity(Mod, Adj);
    WC = WC + pdist(NQ * Ind,'hamming');
end
%% Generate clustering final result
Den = linkage(WC,'average');
if (showden == 1) 
   dendrogram(Den); 
end
FinalResult = cluster(Den, k);