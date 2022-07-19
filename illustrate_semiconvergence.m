% illustrate_semiconvergence.m
%
% Generates Figure 2.1 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% The example illustrates semiconvergence behavior by applying an iterative
% method to an unregularized problem.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May 2021

rng(100)

PSF = psfNSGauss([64, 64], 4, 2, 1.5);
n = 256;

PSF_var = PSF/sum(PSF(:));
optblur.PSF = PSF_var;
[A, b, x, ProbInfo] = PRblur(n, optblur);
[bn,NoiseInfo] = PRnoise(b);

%% Run iterative method to show semiconvergence
opt.x_true = x;
opt.NoStop = 'on';
opt.RegParam = 0;
opt.MaxIter = 300;
iter = 1:300; % return all solutions
[x_lsqr, infogk] = IRhybrid_lsqr(A, bn, iter, opt);

%% select three values of k for illustration
[~,idx] = min(infogk.Enrm);
x_lsqr_min = x_lsqr(:,idx);

idx_small = 10; % too small
x_lsqr_small = x_lsqr(:,idx_small);

idx_large = 200; % too large
x_lsqr_large = x_lsqr(:,idx_large);

%% Get figure of relative errors with reconstructions
c1 = [ 0    0.4470    0.7410];
c2 = [ 0.9290    0.6940    0.1250];
c3 = [ 0.4940    0.1840    0.5560];
c4 = [ 0.4660    0.6740    0.1880];
c5 = [ 0.8500    0.3250    0.0980];
figure,
subplot(2,1,1)
yyaxis left
plot(infogk.Enrm,'LineWidth',3, 'color',c1), hold on
ylabel('rel error norm')

plot(idx_small,infogk.Enrm(idx_small),'*','MarkerSize',12,'LineWidth',2,'color',c2)
plot([idx_small idx_small],[0, .6],'--','color',c2)
plot(idx,infogk.Enrm(idx),'*','MarkerSize',12,'LineWidth',2,'color',c3)
plot([idx idx],[0, .6],'--','color',c3)
plot(idx_large,infogk.Enrm(idx_large),'*','MarkerSize',12,'LineWidth',2,'color',c4)
plot([idx_large idx_large],[0, .6],'--','color',c4)
set(gca,'fontsize',14)

yyaxis right
plot(infogk.Rnrm,'LineWidth',3, 'color',c5), hold on
ylabel('rel residual norm')
xlabel('iteration')

subplot(2,3,4),  imagesc(reshape(x_lsqr_small,n,n)), colormap gray, axis image, axis off
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',c2, 'LineWidth',3);
title('k = 10')
subplot(2,3,5),  imagesc(reshape(x_lsqr_min,n,n)), colormap gray, axis image, axis off
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',c3, 'LineWidth',3);
title(['k = ', num2str(idx)])
subplot(2,3,6),  imagesc(reshape(x_lsqr_large,n,n)), colormap gray, axis image, axis off
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',c4, 'LineWidth',3);
title('k = 200')

