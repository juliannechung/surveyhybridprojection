% generate_tomo.m
%
% Generates Figures 1.2 and 1.3 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%
% The example generates the tomography example.
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May 2021

n=256;
opt = PRtomo('defaults');
opt.angles = 1:179;
[A, b, x, ProbInfo] = PRtomo(n, opt);
bn = PRnoise(b);

figure, imagesc(reshape(x,ProbInfo.xSize)), colormap gray, axis image, axis off
title('true image')
figure, imagesc(reshape(bn,ProbInfo.bSize)), colormap gray, axis image, axis off
title('observations')

%% Get (approximate) inverse solution
opt.NoStop = 'on';
opt.MaxIter = 2000;
opt.RegParam = 0;
[Xinv, infotomo] = IRhybrid_lsqr(A, bn, opt);

figure, imagesc(reshape(Xinv,ProbInfo.xSize)), colormap gray, axis image, axis off
%% Get regularized solution
opt = PRtomo('defaults')
opt.RegParam = 'wgcv';
[Xgk_wgcv, infogk_wgcv] = IRhybrid_lsqr(A, bn, opt);
figure, imagesc(reshape(Xgk_wgcv,ProbInfo.xSize)), colormap gray, axis image, axis off

