clc;
clear;
close all;

addpath('./utils');
% name = 'jaffe';
name = 'Yale';
% name = 'ORL';
% name = 'Isolet';
name = 'usps_random_1000';
% name = 'COIL20';
% name = 'AR';
load (name);
rng('default');  %恢复matlab启动时默认的全局随机流

% fea = NormalizeFea(fea,1);
fea = (mapstd(fea'))'; 

nClusts = length(unique(gnd));%unique除去矩阵中的重复元素
NITER=20;

alpha1=1e-1;
alpha2=1e-1;
lambda=1e0;

beta=7e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters=[1e-5,1e-4,1e-3,1e-2,1e-1,1,10];
% lambdas=[1e-6, 1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1000,1e4,1e5];
filename=['./result/',name,'_',datestr(now,30),'.txt'];
fileID = fopen(filename,'w');
accs=[]; nmis=[];
for i = 1 : length(parameters)
    alpha1 = parameters(i);
    for ii =1 : length(parameters)
        alpha2 = parameters(ii);
        for iii =1 : length(parameters)
            lambda = parameters(iii);
            for iiii =1 : length(parameters)
                beta = parameters(iiii);
            
                fprintf('alpha1 : %f, alpha2 : %f,lambda : %f, beta : %f\n', alpha1,alpha2,lambda,beta);
                fprintf(fileID,'alpha1 : %f, alpha2 : %f, lambda : %f, beta : %f\n', alpha1,alpha2,lambda,beta);
                [Z,W,obj]=SWARG(fea',nClusts,alpha1,alpha2,beta,lambda,NITER);
                
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
                result = ClusteringMeasure(gnd, result_label);
                fprintf(fileID,'all the acc values are :%f\n',result(1));
                fprintf(fileID,'all the nmi values are :%f\n',result(2));
            end
           
        end
    end
end
