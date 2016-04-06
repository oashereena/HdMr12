function [ g1 ] = Evalgx1_4var( xvec,c,nt,p,nv,k)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
g1=zeros(nt,1);
X=c.*ones(nv,1);
for i=1:nt
    X(p)=xvec(i);
g1(i)=gxEval4v(X,k);
end

end

