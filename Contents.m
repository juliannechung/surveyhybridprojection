%% Contents.m
%   This folder contains MATLAB code to accompany the paper:
%
%     "Computational Methods for Large Inverse Problems:
%       A Survey on Hybrid Projection Methods"
%             by Julianne Chung (Virginia Tech) and
%                Silvia Gazzola (University of Bath)
%
%   The DEMO codes require the following packages:
%    (1) IRTools by James Nagy, Silvia Gazzola, and Per Christian Hansen
%              https://github.com/jnagy1/IRtools
% 
%    (2) AIR Tools II package from: https://github.com/jakobsj/AIRToolsII
%
% Chung and Gazzola, May 2022

%% DEMOs
%
%   generate_blur.m     Sets up and runs a 2D image deblurring problem
%                           corresponding to the results in 
%                           Figures 1.1 and 1.3
%
%   generate_tomo.m     Sets up and runs a 2D tomography reconstruction
%                           problem corresponding to the results in 
%                           Figures 1.2 and 1.3
%
%   illustrate_semiconvergence.m     
%                       Generates Figure 2.1 in the paper to illustrate 
%                           semiconvergence phenomenon for image deblurring
%
%   Tikhonov_fixed_lambda.m     
%                       Generates Figure 2.2 in the paper to illustrate 
%                           iterative methods for solving the Tikhonov
%                           problem with a fixed regularization parameter
%
%   Tikhonov_adaptive_lambda.m     
%                       Generates Figure 3.1 in the paper to illustrate 
%                           hybrid projeciton methods where an
%                           iterative method is run with an adaptive
%                           regularization parmaeter selection method.
%
%   compare_subspaces.m
%                       Generates Figure 3.2 in the paper to illustrate 
%                           difference basis images for Arnoldi and
%                           Golub-Kahan bidiagonalization
%  
%  compare_regParamChoice.m
%                       Generates part of Figure 3.3 in the paper to 
%                       compare different regularization parameters for
%                       image deblurring
%
%   estimate_uq_deblur.m
%                       Generates Figure 4.3 in the paper to use GKB and
%                       RSVD to estimate the diagonals and the sum of 
%                       elements of the approximate posterior covariance 
%                       matrix for deblurring
%
%   estimate_uq_tomo.m
%                       Generates Figure 4.3 in the paper to use GKB and
%                       RSVD to estimate the diagonals and the sum of 
%                       elements of the approximate posterior covariance 
%                       matrix for tomography
