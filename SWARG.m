function [Z,W,obj]=SWARG(X,c,alpha1,alpha2,beta,lambda,NITER)
% This code is implemented by Kun Jiang

[d,n] = size(X);
E = ones(n,n);
obj=[];
options = [];
options.NeighborMode = 'KNN';
options.k = 30;
options.WeightMode = 'Binary';
S = constructW(X',options);
S = full(S);
D = diag(sum(S,2));
L = D-S;
[F, ~, ev]=eig1(L, c, 0);

W = 1/d*eye(d);
Z=S;
for iter=1:NITER  
    %update Z
    Z = Z.*((X'*W*W*X+alpha2*S*Z)./(X'*W*W*X*Z+alpha1*E+alpha2*D*Z));
    Z = Z*diag(sqrt(1./(diag(Z'*Z)+eps))); %normalize
    
    %update S
    dist1 = L2_distance_1(F',F');
    dist = dist1+alpha2*L2_distance_1(Z,Z)/lambda;
    for i=1:n
        for j=1:n
            S(i,j)=exp(-dist(i,j)/(2*beta))+eps;
        end
        S(i,:)=S(i,:)./sum(S(i,:));
    end
    LS = (S+S')/2;
    L = diag(sum(LS)) - LS;
    
    %update W
    I= eye(n);
    LZ = (I-Z)*(I-Z');
    T = X*LZ*X';%how to compute T with different graph
    temp1 = 0;
    for i = 1 : d
        temp1 = temp1 + 1/(T(i,i));
    end
    for i = 1 : d
        W(i,i) = 1/(T(i,i) * temp1); 
    end
    
    %update F
    [F, ~, ev]=eig1(L, c, 0);
    
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 10e-11
        lambda = 2*lambda;
    elseif fn2 < 10e-11
        lambda = lambda/2;
    else
        break;
    end
    
    %obj
    tran=0;
    for i1=1:n
        for j1=1:n
            tran=tran+S(i1,j1)*log(S(i1,j1));
        end
    end
    
    obj(iter)=trace((W*X*Z-W*X)*(W*X*Z-W*X)')+alpha1*trace(Z*E)+alpha2*trace(Z'*L*Z)+2*lambda*(trace(F'*L*F)+beta*tran);
%     if iter>2
%         if abs(obj(iter)-obj(iter-1))/obj(iter-1)<1e-8
%             break
%         end
%     end

end

end