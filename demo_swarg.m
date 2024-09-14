clc;
clear;
close all;

addpath('./utils');
% load('usps_random_1000.mat');   % 
load('gaussian_2c_noise.mat');   % 
% load ar_swcan_weights
load MNIST_6996
% fea=X;
% gnd=y;

fea = (mapstd(fea'))'; %mapstd对输入向量进行归一化,使它们具有零均值和单位方差

nClusts = length(unique(gnd));%unique除去矩阵中的重复元素
NITER=30;

alpha1=1e-1;
alpha2=1e-5;%1e-1;

lambda=1e1;
beta=1e-5;%beta=7e-3;

tic
[Z,W,obj]=SWARG(fea',nClusts,alpha1,alpha2,beta,lambda,NITER);
toc

% plot(obj);

%ncut
addpath('Ncut_9');
A = Z;
A = (A+A')/2;  
[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(A,nClusts);
result_label = zeros(size(fea,1),1);%vec2ind
for j = 1:nClusts
    id = find(NcutDiscrete(:,j));
    result_label(id) = j;
end
result_n = ClusteringMeasure(gnd, result_label)

figure;plot(W);


