function [partition, m, K, partition_bin] = pycrprnd(theta, sigma, n)

%
% PYCRPRND samples a partition from the Pitman-Yor Chinese restaurant process
%   [partition, m, K, partition_bin] = pycrprnd(alpha, sigma, n)
%
%--------------------------------------------------------------------------
% INPUTS
%   - theta:    scale parameter of the PY CRP
%   - sigma:    discount parameter of the PY CRP
%   - n:        number of objects in the partition
%
% OUTPUTS
%   - partition: Vector of length n. partition(i) is the cluster membership
%               of object i. Clusters are numbered by order of appearance
%   - m:        Vector of length n. m(j) is the size of cluster j.
%               m(j)=0 for j>K
%   - K:        Integer. Number of clusters
%   - partition_bin: locical matrix of size n*n. partition(i,j)=true if object
%               i is in cluster j, false otherwise
%--------------------------------------------------------------------------
% EXAMPLE
% theta = 3; sigma = .5; n= 100;
% [partition, m, K, partition_bin] = pycrprnd(theta, sigma, n);
%--------------------------------------------------------------------------

% Check parameters
if sigma<0 || sigma>=1
    error('Parameter sigma must be in [0,1)');
end
if theta<= -sigma
    error('Parameter theta must be in [-sigma,Inf)');
end

m = zeros(n, 1);
partition = zeros(n, 1);

% Initialization
partition(1) = 1;
m(1) = 1;
K = 1;
% Iterations
for i=2:n
    % Compute the probability of joining an existing cluster or a new one
    proba = [m(1:K)-sigma; theta + sigma*K]/(theta+i-1);
    % Sample from a discrete distribution w.p. proba
    u = rand;
    partition(i) = find(u<=cumsum(proba), 1);
    % Increment the size of the cluster
    m(partition(i)) = m(partition(i)) + 1;
    % Increment the number of clusters if new
    K = K + isequal(partition(i), K+1);    
end
partition_bin = sparse((1:n)', partition, ones(n,1)); % partition_bin(i,j) =1 if item i is in cluster j
m = m(1:K);