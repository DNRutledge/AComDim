%% Load data
% load('Lignin_40.mat');
% IN
% F:\Matlab\toolbox\IAQA

%%
X=Lignine_Starch_smooth_Balanced;
[nR,nC]=size(X);

[nR,nC]=size(X);
SampleNums=[1:nR]';
VarNums=[1:nC];

%% Do A_ComDim to get back the ANOVA matrices
Opts.Plots=2;
Opts.Dims=6;
Opts.ndim=6;

Opts.decomp='classical';
% Opts.decomp='glm';
Opts.interactions=3;

r=3;
c=2;

Opts.Samples=0;
lbl_var=[];
% lbl_var=VarNums;

%%
Opts.MVA = 'PCA';

[A_ComDimIn, A_ComDimOut, A_ComDimResult, A_Collection]=A_ComDim_2020_LMS(X, lbl_var, Grps, Opts);


%%
Opts.MVA = 'ICA';

[A_ComDimIn, A_ComDimOut, A_ComDimResult, A_Collection]=A_ComDim_2020_LMS(X, lbl_var,Grps, Opts);

