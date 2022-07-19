% estimate_uq_tomo.m
%
% Generates Figure 4.3 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% For the tomography example, this script uses GKB to estimate the 
%   diagonals and the sum of elements of the approximate posterior 
%   covariance matrix, as a function of GKB iterations.
%
%   Comparisons to approximations obtained using low-rank RSVD
%       are made.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May, 2022

%% Tomography example
n = 256;
opt = PRtomo('defaults');
opt.angles = 1:179;
[A, b, x, ProbInfo] = PRtomo(n, opt);
[bn,noiseInfo] = PRnoise(b);
sigma = std(noiseInfo.noise);

figure, imagesc(reshape(x,ProbInfo.xSize)), colormap gray, axis image, axis off
title('true image')
figure, imagesc(reshape(bn,ProbInfo.bSize)), colormap gray, axis image, axis off
title('observations')

%% get regularized solution
opt = PRtomo('defaults');
opt.RegParam = 'wgcv';
[~, infogk_wgcv_auto] = IRhybrid_lsqr(A, bn, opt);

%% For IRhybrid_lsqr to run to MaxIter
opt = PRtomo('defaults');
opt.RegParam = 'wgcv';
opt.Reorth = 'on';
opt.DecompOut = 'on';
opt.NoStop = 'on';
opt.MaxIter = 100;
[Xgk_wgcv, infogk_wgcv] = IRhybrid_lsqr(A, bn, opt);
figure, imagesc(reshape(Xgk_wgcv,ProbInfo.xSize)), colormap gray, axis image, axis off
title('regularized solution')

%% UQ - GKB estimates for various iterations
Bk = infogk_wgcv.B;
Vk = infogk_wgcv.V;
lambda = infogk_wgcv.RegP(end)^2;

k = size(Bk,2);
d = zeros(n^2,k);
dsum = zeros(1,k);
for i = 1:k
   [d(:,i), dsum(i)] = diaginv(Bk(1:i+1,1:i), Vk(:,1:i), lambda, sigma);
end

%% UQ - RSVD estimates for various ranks
d_rsvd = zeros(n^2,k);
dsum_rsvd = zeros(1,k);
for i = 1:k
    [d_rsvd(:,i),dsum_rsvd(i)] = diaginv_rsvd(A, lambda, sigma, i);
end

%% create figures of variance estimates (diagonal entries)
iter = infogk_wgcv_auto.its;

bottom = min(min(d(:,iter)),min(min(d_rsvd(:,iter))));
top = max(max(d(:,iter)),max(max(d_rsvd(:,iter))));

figure, 
tlt = tiledlayout(1,2);
nexttile
imagesc(reshape(d(:,iter),ProbInfo.xSize)), 
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',colors(1), 'LineWidth',3);

colormap gray, axis image, axis off
caxis manual
caxis([bottom top]);
title('GKB'), 
set(gca,'fontsize',22)

nexttile
imagesc(reshape(d_rsvd(:,iter),ProbInfo.xSize)),
hold on, plot([0 0 256  256 0 ],[0 256 256 0 0],'color',colors(2), 'LineWidth',3);
colormap gray, axis image, axis off
caxis manual
caxis([bottom top]);
title('RSVD'), 
set(gca,'fontsize',22)

tlt.Padding = "tight";
tlt.TileSpacing = "tight";

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Ruler.Exponent = -3;
xlength = 600;
ylength = 400;
set(gcf,'Position',[100 100 xlength ylength])
% saveas(gcf,'figs/cov_tomo.eps','epsc')

%% plot sum of covariance-variance estimates
figure, semilogy(dsum,'-', 'LineWidth',2, 'color',colors(1)), hold on,
semilogy(dsum_rsvd,'--', 'LineWidth',2, 'color',colors(2))
semilogy([iter iter],[10^-1 10^3],':k')
legend('GKB','RSVD','')

ylabel('sum of variance-covariance')
xlabel('iteration')
axis([0 opt.MaxIter 10^-1 10^3])
set(gca,'fontsize',22)

xlength = 600;
ylength = 400;
set(gcf,'Position',[100 100 xlength ylength])
% saveas(gcf,'figs/sumcov_tomo.eps','epsc')
