% compare_regParamChoice

% Generates part of Figure 3.3 in the paper 
%   "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"

% This script compares image deblurring for
%   different regularization parameters 
%
% Silvia Gazzola, University of Bath
% Julianne Chung, Virginia Tech
% May, 2021

n = 256; % size of the problem

%% image deblurring example
PSF = psfNSGauss([64, 64], 4, 2, 1.5);
% rescaling the PSF values -- meaningful if we think of the physical
% interpretation of PSF
PSF_var = PSF/sum(PSF(:));
optblur.PSF = PSF_var;
% generating a blurring problem with this PSF
[A, b, x] = PRblur(n, optblur);

% setting options for solvers
MaxIter = 100;

opt.x_true = x;
opt.NoStop = 'on';
opt.MaxIter = MaxIter;

opt1 = opt;
opt1.RegParam = 'optimal';
opt2 = opt;
opt2.RegParam = 'discrepit';
opt2.NoiseLevel = 1e-2;
opt3 = opt;
opt3.RegParam = 'wgcv';
opt3.stopGCV = 'GCVvalues';
opt3.GCVflatTol = 10^-16;

%% creating histogram plots
tries = 100;

% allocating memory
RegPopt = zeros(MaxIter, tries);
Enrmopt = zeros(MaxIter, tries);
BestIt = zeros(tries,1);
BestRegP = zeros(tries,1);
BestEnrm = zeros(tries,1);

RegPdp = zeros(MaxIter, tries);
Enrmdp = zeros(MaxIter, tries);
StopItdp = zeros(tries,1);
StopRegPdp = zeros(tries,1);
StopEnrmdp = zeros(tries,1);
BestItdp = zeros(tries,1);
BestRegPdp = zeros(tries,1);
BestEnrmdp = zeros(tries,1);

RegPwgcv = zeros(MaxIter, tries);
Enrmwgcv = zeros(MaxIter, tries);
StopItwgcv = zeros(tries,1);
StopRegPwgcv = zeros(tries,1);
StopEnrmwgcv = zeros(tries,1);
BestItwgcv = zeros(tries,1);
BestRegPwgcv = zeros(tries,1);
BestEnrmwgcv = zeros(tries,1);

for i = 1:tries
    disp(i)
    bn = PRnoise(b);
    % dp
    [~, infogk2temp] = IRhybrid_lsqr(A, bn, opt2);
    % wgcv
    [~, infogk3temp] = IRhybrid_lsqr(A, bn, opt3);
    %
    RegPdp(:,i) = infogk2temp.RegP;
    RegPwgcv(:,i) = infogk3temp.RegP;
    Enrmdp(:,i) = infogk2temp.Enrm;      
    Enrmwgcv(:,i) = infogk3temp.Enrm;   
    %
    StopItdp(i) = infogk2temp.StopReg.It;
    StopItwgcv(i) = infogk3temp.StopReg.It;
    StopRegPdp(i) = infogk2temp.StopReg.RegP;
    StopRegPwgcv(i) = infogk3temp.StopReg.RegP;
    StopEnrmdp(i) = infogk2temp.StopReg.Enrm;
    StopEnrmwgcv(i) = infogk3temp.StopReg.Enrm;
    %
    BestItdp(i) = infogk2temp.BestReg.It;
    BestRegPdp(i) = infogk2temp.BestReg.RegP;
    BestEnrmdp(i) = infogk2temp.BestReg.Enrm;
    %
    BestItwgcv(i) = infogk3temp.BestReg.It;
    BestRegPwgcv(i) = infogk3temp.BestReg.RegP;
    BestEnrmwgcv(i) = infogk3temp.BestReg.Enrm;   
end

plot2_dp = [StopItdp, StopRegPdp];
plot2_wgcv = [StopItwgcv, StopRegPwgcv];

figure
hist3(plot2_dp,'Ctrs',{12:1:25 [0.029:0.0001:0.034, 0.05, 0.071:0.0001:0.073]}, 'FaceColor', 'r')
hold on
hist3(plot2_wgcv,'Ctrs',{12:1:25 [0.029:0.0001:0.034, 0.05, 0.071:0.0001:0.073]}, 'FaceColor', 'y')
set(gca,'FontSize',20)
xticks(12:2:25)
xlabel('k')
ylabel('\lambda_k')
legend('DP', 'wGCV')

%% creating surface plots
nl = 100;
lambdasample = logspace(-6,1,nl);

EnrmSurf = zeros(MaxIter, nl);
for i = 1:nl
    disp(i)
    opt.RegParam = lambdasample(i);
    [Xgktemp, infogktemp] = IRhybrid_lsqr(A, bn, opt);
    EnrmSurf(:,i) = infogktemp.Enrm;
end

[Xgkdp, infogkdp] = IRhybrid_lsqr(A, bn, opt2);
[Xgkwgcv, infogkwgcv] = IRhybrid_lsqr(A, bn, opt3);

% % error surface
[XX,YY] = meshgrid(lambdasample,1:MaxIter);
figure, mesh(log10(XX), YY, log10(EnrmSurf))
view([-58.0 61])
xlabel('log_{10}(\lambda_k)')
ylabel('k')
zlabel('log_{10}(RRE)')
set(gca,'FontSize',22)

% % aerial view
[XX,YY] = meshgrid(lambdasample,1:MaxIter);
figure, contourf(log10(XX),YY,log10(EnrmSurf))
hold on
plot(log10(infogkdp.RegP), 1:MaxIter, '*', 'MarkerSize', 10, 'LineWidth', 2, 'color', 'r')
plot(log10(infogkwgcv.RegP), 1:MaxIter, '+', 'MarkerSize', 10, 'LineWidth', 2, 'color', 'y')
colorbar
xlabel('log_{10}(\lambda_k)')
ylabel('k')
legend('log_{10}(error surface)', 'discrepancy', 'wgcv')
legend('Location','northwest')
set(gca,'FontSize',22)



