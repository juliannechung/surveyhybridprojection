function [d, sm] = diaginv(Bk, Vk, lambda, sigma)
%   [d, sm] = diaginv(Bk, Vk, lambda, sigma)
%
% This function computes the diagonals and the sum of elements of 
%   an approximate posterior covariance matrix
%           \hat{Gamma_post} \approx Gamma_post
%   where
%   \hat{Gamma_post} = sigma^2 (\lambda I + V_k B_k^T B_k V_k^T)^{-1}
%
%   Input:
%       Bk - (k+1)xk bidiagonal matrix from Golub-Kahan bidiagonalization
%       Vk - nxk orthonormal basis vectors from Golub-Kahan bidiagonalization
%       lambda - regularization parameter
%       sigma - noise standard deviation
%
%   Output:
%       d - vector of diagonal entries
%       sm - sum of all entries in \hat{Gamma_post}
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May, 2022

% (1) estimate the diagonals
[~,s,W] = svd(Bk,0);
s = diag(s);
d12 = s./sqrt(s.^2 + lambda); % diagonal elements of Delta_k^.5 (without lambda)
invL = (Vk*W)*diag(d12);

d =  ones(size(invL,1),1) - diaglowrank(invL);
d = d./lambda; % include lambda
d = sigma^2 * d; % vector of diagonal entries

% (2) estimate the sum of variance-covariance entries
n = size(Vk,1);
w = W'*(Vk'*ones(n,1));

sw = w.*d12;
sm =  n - sw'*sw;
sm = sm./lambda;
sm = sigma^2 * sm;

end

function d = diaglowrank(V)
% computes the diagonal elements of a low rank matrix V V'

n = size(V,1);
d = zeros(n,1);
for i = 1:n
    d(i) = V(i,:)*V(i,:)';
end

end
