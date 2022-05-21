%WISED SPECTRAIL Generate individual Spectral Clustering.
%
%   Input  : adj          : adjancency matrix

%            k            : desired number of nodes in groups [n1, n2, ..]
%            sigma        : for converting the sparse distance matrix to a sparse similarity matrix
   
%
%   Output : cluster_labels: partitional result of clustering,
%            modules       : modular result of clustering.
%
% Written by: Muhammad Yousefnezhad,
%             iBRAIN group, Department of Computer Science & Technology,
%             Nanjing University of Aeronautics & Astronautics  
function [cluster_labels, modules] = wisedSpectral(adj, k, sigma)
%% Convert the sparse distance matrix to a sparse similarity matrix,
% where S = exp^(-(A^2 / 2*sigma^2)).
% Note: This step can be ignored if A is sparse similarity matrix.
disp('Converting distance matrix to similarity matrix...');
A = adj;
n = size(A, 1);
if (sigma == 0) % Selftuning spectral clustering
  % Find the count of nonzero for each column
  disp('Selftuning spectral clustering...');
  col_count = sum(A~=0, 1)';
  col_sum = sum(A, 1)';
  col_mean = col_sum ./ col_count;
  [x y val] = find(A);
  A = sparse(x, y, -val.*val./col_mean(x)./col_mean(y)./2);
  clear col_count col_sum col_mean x y val;
else % Fixed-sigma spectral clustering
  disp('Fixed-sigma spectral clustering...');
  A = A.*A;
  A = -A/(2*sigma*sigma);
end
%% Do exp function sequentially because of memory limitation
num = 2000;
num_iter = ceil(n/num);
S = sparse([]);
for i = 1:num_iter
  start_index = 1 + (i-1)*num;
  end_index = min(i*num, n);
  S1 = spfun(@exp, A(:,start_index:end_index)); % sparse exponential func
  S = [S S1];
  clear S1;
end
clear A;
%% Do laplacian, L = D^(-1/2) * S * D^(-1/2)
disp('Doing Laplacian...');
D = sum(S, 2) + (1e-10);
D = sqrt(1./D); % D^(-1/2)
D = spdiags(D, 0, n, n);
L = D * S * D;
%% Do eigendecomposition 
% if L =
%   D^(-1/2) * S * D(-1/2)    : set 'LM' (Largest Magnitude), or
%   I - D^(-1/2) * S * D(-1/2): set 'SM' (Smallest Magnitude).
%
disp('Performing eigendecomposition...');
OPTS.disp = 0;
[Vn, ~] = eigs(L, k, 'LM', OPTS);
%% Normalize each row to be of unit length
disp('Performing kmeans...');
sq_sum = sqrt(sum(Vn .* Vn, 2)) + 1e-20;
U = Vn ./ repmat(sq_sum, 1, k);
%% Do Cluserting
cluster_labels = k_means(U, [], k, 3);
%% Do Modularity
% find the Fiedler vector: eigenvector corresponding 
% to the second smallest eigenvalue of the Laplacian matrix
S = full(S);
[V,D2] = eig(full(L));
[ds,Y] = sort(diag(D2));
fv = V(:,Y(2));
[~,I]=sort(fv);
k = [0 k];
for kk=1:length(k)
modules{kk}=[];
for x=1:k(kk); modules{kk} = [modules{kk} I(x+k(kk-1))]; end
end
modules = modules(2:length(modules));
end