% generate_blur.m
%
% Generates Figures 1.1 and 1.3 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% The example generates the image deblurring example.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May 2021
rng(0)

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
figure, imagesc(reshape(x,n,n)), colormap gray, axis image, axis off
title('true image')
figure, imagesc(reshape(bn,n,n)), colormap gray, axis image, axis off
title('observed image')
figure, imagesc(PSFshow), colormap gray, axis image, axis off
title('PSF')

%% approximate the blurring matrix and compute naive and TSVD solutions
A1 = Kall.a{1};
B1 = Kall.b{1};
B = B1*(reshape(x, n, n)*A1');
B = PRnoise(B);
figure, imagesc(B), colormap gray, axis image, axis off
Xnaive = B1\(B/A1');
figure, imagesc(Xnaive), colormap gray, axis image, axis off
[UA1, SA1, VA1] = svd(A1);
[UB1, SB1, VB1] = svd(B1);
trunc = 35;
Xreg = (VB1(:,1:trunc)*diag(1./diag(SA1(1:trunc,1:trunc)))*UB1(:,1:trunc)')*(B*(UA1(:,1:trunc)*diag(1./diag(SA1(1:trunc,1:trunc)))*VA1(:,1:trunc)'));

figure, imagesc(Xreg), colormap gray, axis image, axis off
% % try to run some solvers
opt.x_true = x;
opt.NoStop = 'on';
[Xgk1, infogk1] = IRhybrid_lsqr(A, bn, opt);

figure, imagesc(reshape(Xnaive,ProbInfo.xSize)), colormap gray, axis image, axis off
title('inverse solution')
figure, imagesc(reshape(Xreg,ProbInfo.xSize)), colormap gray, axis image, axis off
title('regularized solution')