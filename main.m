clc;
clear;
close all;

addpath('./utils');
% load('AR10P.mat');   % 1e-2 1e1 1e0 1e-3
load('Yale.mat');   % 
% load('COIL20.mat');   % 


% fea = NormalizeFea(fea,1);
fea = (mapstd(fea'))';

nClusts = length(unique(gnd));%unique除去矩阵中的重复元素
NITER=50;

alpha1=0.000010;
alpha2=0.001000;
lambda=0.000010;
beta=0.001000;
%for no weights

alpha1=1e-2;
alpha2=1e-1;
lambda=1e-1;
beta=5e-3;

tic
[Z,W,obj]=SWARG(fea',nClusts,alpha1,alpha2,beta,lambda,NITER);
toc


addpath('Ncut_9');
A = Z;
A = A - diag(diag(A));
A = (A+A')/2;  

[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(A,nClusts);
result_label = zeros(size(fea,1),1);%vec2ind
for j = 1:nClusts
    id = find(NcutDiscrete(:,j));
    result_label(id) = j;
end
result = ClusteringMeasure(gnd, result_label)

plot(obj);

