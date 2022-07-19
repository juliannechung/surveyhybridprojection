function [d, sest] = diaginv_rsvd(A, lambda, sigma, l)
%   [d, sest] = diaginv_rsvd(A, lambda, sigma, l)
%
% This function computes the diagonals and the sum of elements of 
%   an approximate posterior covariance matrix
%           \hat{Gamma_post} \approx Gamma_post
%   where \hat{Gamma_post} is constructed using a low-rank approximation
%   based on the randomized SVD (RSVD)
%
%   Input:
%       A - matrix or function handle
%       lambda - regularization parameter
%       sigma - noise standard deviation
%       l - dimension of RSVD approximation
%
%   Output:
%       d - vector of diagonal entries
%       sm - sum of all entries in \hat{Gamma_post}
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May, 2022

[~,s,W] = getRSVD(A,l,1); % RSVD approximation
s = diag(s);

% (1) diagonal estimation
d12 = s./sqrt(s.^2 + lambda);
invL = W*diag(d12);

d =  ones(size(invL,1),1) - diaglowrank(invL);
d = d./lambda; % final scaling
d = sigma^2 * d;

% (2) sum of all entries
n = size(A,2);
w = W'*ones(size(W,1),1);
sw = w .* d12;
sest =  n - sw'*sw;
sest = sest./lambda; % final scaling
sest = sigma^2 * sest;


end

function [U,S,V] = getRSVD(A, l, flag)
% Generates RSVD of rank l
%   Inputs: A - input matrix
%           l - dimension of approximation
%           flag - 0 left multiplication with random matrix
%                  1 right multiplication with random matrix
%
%   Outputs: U,S,V - components of the RSVD
%
if nargin < 3
    flag = 0;
end

switch flag
    case 0        
        m = size(A,1);
        W = randn(l, m);
        for i = 1:size(W,1)
%         Y = W*A;
            Y(:,i) = A'*W(i,:)';
        end
        [Q,~] = qr(Y,0);
        for i = 1:size(Q,2)
            B(:,i) = A*Q(:,i);
        end
        
        [U, S, H] = svd(B, 'econ');
        V = Q*H;
        
    case 1
        % Good for approximating V
        n = size(A,2);
        W = randn(n, l);
        for i = 1:l
            Y(:,i) = A*W(:,i);
        end
        [Q,~] = qr(Y,0);
        %         B = Q'*A;
        for i = 1:size(Q,2)
           B(:,i) = A'*Q(:,i);
        end
        [O, S, V] = svd(B', 'econ');
        U = Q*O;
end

end

function d = diaglowrank(V)
% computes the diagonal elements of a low rank matrix V V'
n = size(V,1);
d = zeros(n,1);
for i = 1:n
    d(i) = V(i,:)*V(i,:)';
end

end


