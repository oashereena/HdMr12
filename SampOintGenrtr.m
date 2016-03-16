function[xvec]=SampOintGenrtr(xmax,xmin,nt)
mux=(xmax+xmin)/2;
sx=(xmax-xmin)/(nt-1);
xvec=zeros(5,1);
for i=1:nt
xvec(i)=mux+sx*(2*i-nt-1)/2;
end