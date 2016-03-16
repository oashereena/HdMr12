function[PHI]=LagrangeInterpolator(xsp,nt,x)
PHI = sym('phi',[nt,1]);
PHI(1:end,1:end)=zeros(nt,1);
for i=1:nt
xvdummy=[xsp(1:i-1);xsp(i+1:end)];
PHI(i)=symfun((prod(x-xvdummy))/(prod(xsp(i)-xvdummy)),x);
end