function [Z,W,obj]=SWARG_RC_PG(X,c,alpha1,NITER)
% This code is implemented by Kun Jiang

[d,n] = size(X);
E = ones(n,n);
obj=[];
options = [];
options.NeighborMode = 'KNN';
options.k = 30;
options.WeightMode = 'Binary';
S = constructW(X',options);
Z = full(S);

W = 1/d*eye(d);
for iter=1:NITER  
    %update Z
    Z = Z.*((X'*W*W*X)./(X'*W*W*X*Z+alpha1*E));
    Z = Z*diag(sqrt(1./(diag(Z'*Z)+eps))); %normalize

    
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
    
    obj(iter)=trace((W*X*Z-W*X)*(W*X*Z-W*X)')+alpha1*trace(Z*E);
    if iter>2
        if abs(obj(iter)-obj(iter-1))/obj(iter-1)<1e-8
            break
        end
    end

end

end