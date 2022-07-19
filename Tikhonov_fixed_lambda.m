% Tikhonov_fixed_lambda.m
%
% Generates Figure 2.2 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% The example illustrates applying an iterative method to a Tikhonov
% problem, with fixed regularization.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May 2021

rng(100)

PSF = psfNSGauss([64, 64], 4, 2, 1.5);
PSFshow = PSF(16+1:16+32, 16+1:16+32);
PSFshow = imresize(PSFshow,2);
n = 256;

PSF_var = PSF/sum(PSF(:));
optblur.PSF = PSF_var;

[A, b, x, ProbInfo] = PRblur(n, optblur);
[bn,NoiseInfo] = PRnoise(b);

%% Run iterative method with different fixed lambdas
opt.x_true = x;
opt.NoStop = 'on';
opt.MaxIter = 300;

opt.RegParam = 0;
[x_lsqr, infogk] = IRhybrid_lsqr(A, bn, opt);

opt.RegParam = 5e-3;
[Xgk1, infogk1] = IRhybrid_lsqr(A, bn, opt);

opt.RegParam = 2e-2;
[Xgk2, infogk2] = IRhybrid_lsqr(A, bn, opt);

opt.RegParam = 5e-1;
[Xgk3, infogk3] = IRhybrid_lsqr(A, bn, opt);

%% Get figure of relative errors with reconstructions
c1 = [0    0.4470    0.7410];
c2 = [0.9290    0.6940    0.1250]; 
c3 = [0.4940    0.1840    0.5560];
c4 = [0.4660    0.6740    0.1880];
figure,
subplot(2,1,1)
plot(infogk.Enrm,'-', 'LineWidth',2, 'color',c1), hold on
plot(infogk1.Enrm,'--', 'LineWidth',2, 'color',c2)
plot(infogk2.Enrm,'-.', 'LineWidth',2, 'color',c3)
plot(infogk3.Enrm,':', 'LineWidth',2, 'color',c4)
ylabel('rel error norm')
xlabel('iteration')
axis([0 300 .14 .4])
set(gca,'fontsize',14)

plot(opt.MaxIter,infogk1.Enrm(opt.MaxIter),'*','MarkerSize',12,'LineWidth',2,'color',c2)
plot(opt.MaxIter,infogk2.Enrm(opt.MaxIter),'*','MarkerSize',12,'LineWidth',2,'color',c3)
plot(opt.MaxIter,infogk3.Enrm(opt.MaxIter),'*','MarkerSize',12,'LineWidth',2,'color',c4)
legend('\lambda=0', '\lambda=5e-3', '\lambda=2e-2', '\lambda=5e-1','Location','nw')

% undersmoothed/noisy
subplot(2,3,4),  imagesc(reshape(Xgk1,n,n)), colormap gray, axis image, axis off
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',c2, 'LineWidth',3);
title('\lambda = 5e-3')
% ok
subplot(2,3,5),  imagesc(reshape(Xgk2,n,n)), colormap gray, axis image, axis off
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',c3, 'LineWidth',3);
title('\lambda = 2e-2')
% oversmoothed
subplot(2,3,6),  imagesc(reshape(Xgk3,n,n)), colormap gray, axis image, axis off
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',c4, 'LineWidth',3);
title('\lambda = 5e-1')


