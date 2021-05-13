% Tikhonov_adaptive_lambda.m
%
% Generates Figure 2.3 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% The example illustrates applying an iterative method to a Tikhonov
% problem, with adaptive regularization parameter.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% April 2021

rng(100)

PSF = psfNSGauss([64, 64], 4, 2, 1.5);
PSFshow = PSF(16+1:16+32, 16+1:16+32);
PSFshow = imresize(PSFshow,2);
n = 256;

PSF_var = PSF/sum(PSF(:));
optblur.PSF = PSF_var;

[A, b, x, ProbInfo] = PRblur(n, optblur);
[bn,NoiseInfo] = PRnoise(b);

%% Trying some automatic parameter choice strategies
opt.x_true = x;
opt.NoStop = 'on';
opt.MaxIter = 300;

opt.RegParam = 2e-2;
[x_best, info_best] = IRhybrid_lsqr(A, bn, opt);

opt.RegParam = 'discrepit';
opt.NoiseLevel = 1e-2;
[Xgk_dp, infogk_dp] = IRhybrid_lsqr(A, bn, opt);

opt.RegParam = 'wgcv';
[Xgk_wgcv, infogk_wgcv] = IRhybrid_lsqr(A, bn, opt);

opt.RegParam = 'optimal';
[Xgk_opt, infogk_opt] = IRhybrid_lsqr(A, bn, opt);

%% Get figure of relative errors with reconstructions

c1 = [0.3010    0.7450    0.9330];
c2 = [0.8500    0.3250    0.0980]; 
c3 = [0.4940    0.1840    0.5560];

figure,
subplot(1,2,1)
plot(infogk_dp.Enrm,'-', 'LineWidth',2, 'color',c1), hold on
plot(infogk_wgcv.Enrm,'--', 'LineWidth',2, 'color',c2)
plot(infogk_opt.Enrm,'-.', 'LineWidth',2, 'color',c3), 
ylabel('rel error norm')
xlabel('iteration')
set(gca,'fontsize',14)
legend('adaptive DP','adaptive wGCV','adaptive opt','Location','ne')

subplot(1,2,2)
semilogy(infogk_dp.RegP,'-', 'LineWidth',2, 'color',c1), hold on
semilogy(infogk_wgcv.RegP,'--', 'LineWidth',2, 'color',c2), hold on
semilogy(infogk_opt.RegP,'-.', 'LineWidth',2, 'color',c3)
semilogy([1 300],[5e-3 5e-3],':k')
semilogy([1 300],[2e-2 2e-2],':k')
semilogy([1 300],[5e-1 5e-1],':k')
yticks([1e-5 5e-3 2e-2 5e-1 1.1])
ylabel('\lambda_k')
yticklabels({'1e-5','5e-3','2e-2','5e-1','1.1'})
xlabel('iteration')
set(gca,'fontsize',14)

