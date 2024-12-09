function W=gpi(A,B,s)
if nargin<3
    s=1;
end
[m,k]=size(B);
if m<k 
    disp('Warning: error input!!!');
    W=null(m,k);
    return;
end
A=max(A,A');

if s==0
    alpha=abs(eigs(A,1));
else if s==1
    ww=rand(m,1);
    for i=1:10
        m1=A*ww;
        q=m1./norm(m1,2);
        ww=q;
    end
    alpha=abs(ww'*A*ww);
    else disp('Warning: error input!!!');
         W=null(m,k);
         return;
    end
end
        
err1=1;t=1;
W=orth(rand(m,k));
A_til=alpha.*eye(m)-A;
while t<5
    M=2*A_til*W+2*B;
    [U,~,V]=svd(M,'econ');
    W=U*V';
    obj(t)=trace(W'*A*W-2.*W'*B);
    if t>=2
        err1=abs(obj(t-1)-obj(t));
    end
        t=t+1;
end