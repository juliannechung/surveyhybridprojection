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
% opt.angles = 1:0.8:179;
opt.angles = 1:179;
[A, b, x, ProbInfo] = PRtomo(n, opt);
bn = PRnoise(b);

figure, imagesc(reshape(x,ProbInfo.xSize)), colormap gray, axis image, axis off
title('true image')
figure, imagesc(reshape(bn,ProbInfo.bSize)), colormap gray, axis image, axis off
title('observations')

%% Get inverse solution
warning('Getting the inverse solution can take a long time')
[Xinv] = lsqr(A,bn,[],64798);

figure, imagesc(reshape(Xinv,ProbInfo.xSize)), colormap gray, axis image, axis off
title('inverse solution')

%% get regularized solution
opt = PRtomo('defaults')
opt.RegParam = 'wgcv';
[Xgk_wgcv, infogk_wgcv] = IRhybrid_lsqr(A, bn, opt);
figure, imagesc(reshape(Xgk_wgcv,ProbInfo.xSize)), colormap gray, axis image, axis off
title('regularized solution')
