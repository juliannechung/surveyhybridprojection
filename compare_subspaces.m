% compare_subspaces.m
%
% Generates Figure 3.2 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% This script compares Arnoldi and GK subspaces for the image deblurring
% example and provides GK subspaces for the tomography example.  Subspace
% vectors are reshaped into images.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May 2021

%% image deblurring example
PSF = psfNSGauss([64, 64], 4, 2, 1.5);
% rescaling the PSF values -- meaningful if we think of the physical
% interpretation of PSF
PSF = PSF/sum(PSF(:));
% generating a blurring problem with this PSF
optblur.PSF = PSF;
PSFshow = PSF(16+1:16+32, 16+1:16+32);
PSFshow = imresize(PSFshow,2);
n = 256;
[A, b, x, ProbInfo] = PRblur(n, optblur);
Kall = kronApprox(A);
bn = PRnoise(b);

% Run hybrid projection methods
opt.x_true = x;
opt.NoStop = 'on';
opt.DecompOut = 'on';

[Xgm1, infogm1] = IRhybrid_gmres(A, bn, opt);

[Xgk1, infogk1] = IRhybrid_lsqr(A, bn, opt);

% Get subspace vectors
V_A = infogm1.V;
V_G = infogk1.V;

%% Tomography example
opt = PRtomo('defaults');
opt.angles = 1:179;
[A, b, x, ProbInfo] = PRtomo(n, opt);
bn = PRnoise(b);

opt.x_true = x;
opt.NoStop = 'on';
opt.DecompOut = 'on';
[Xtomo, infotomo] = IRhybrid_lsqr(A, bn, opt);
V_tomo = infotomo.V;

%% Compare solution subspace vectors as images
idx = [1 2 4 10];
figure
for i = 1:length(idx)
    subplot(3, length(idx),i), imshow(reshape(V_A(:,idx(i)),n,n),[]), colormap gray, axis image, axis off
    title(sprintf('k = %d',idx(i))), set(gca,'fontsize',14)
    subplot(3, length(idx),length(idx)+i), imshow(reshape(V_G(:,idx(i)),n,n),[])
    subplot(3, length(idx),2*length(idx)+i), imshow(reshape(V_tomo(:,idx(i)),n,n),[])
end
subplot(3,length(idx),1), ylabel('deblur (Ar)'),set(gca,'fontsize',14)
subplot(3,length(idx),length(idx)+1), ylabel('deblur (GK)'),set(gca,'fontsize',14)
subplot(3,length(idx),2*length(idx)+1), ylabel('tomo (GK)'),set(gca,'fontsize',14)

return
