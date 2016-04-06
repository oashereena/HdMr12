% Cantilever beam with tip point load P, flexural rigidity EI, length l
% y=(Pl^3)/(3EI) modeled as x=k(theta1)(theta2)
% X as tip deflection as measurements, in mm X1,X2, ... observations
% theta1 denote as th1 as P ~ N(10,1) kN, Bounds (9,11)kN
% theta2 denote as th2 as (1/EI) ~ N(0.0007,0.0002) kNm2 , Bounds (0.0006,0.0008)kNm2
% k = (1000/3) =333.3333, with l= 1m, and with adjustments in units made for y to be in mm , included
% symbolic integration x for theta1 and y for theta2
% nrfTHX numerator term of f theta conditioned on X ;drfTHx for denominator
% mDx marginal distribution of theta1 ; mDy marginal distribution of theta2
% xhat estimate of theta1; yhat estimate of theta2
k=1000/3;muth1=10.0;muth2=1500;sigth1=1;sigth2=500;th1min=9;th1max=11;th2min=1000;th2max=2000;X1=2.9202;
syms x y
muX=k*muth1*vpa((int((1/x)*normpdf(x,muth2,sigth2),x,th2min,th2max)),5);
c1=double(vpa(muX/k,5));
Ig=@(x,y) ((x/y)-c1).^2*exp(-0.5*((x-muth1).^2+((y-muth2)./sigth2).^2));
sig2X=k^2*(1/(2*pi*sigth1*sigth2))*vpa(int(int(Ig(x,y),x,th1min,th1max),y,th2min,th2max),5);
sX=sqrt(sig2X);
mux=muth1;muy=muth2;sx=sigth1;sy=sigth2;xmin=th1min;xmax=th1max;ymin=th2min;ymax=th2max;
syms x y
nrfTHX=exp(-0.5*(((x-mux)/sx)^2+((y-muy)/sy)^2+((X1-k*x*y^-1)/sX)^2));
nrfTHXh=@(x,y) exp(-0.5.*(((x-double(mux))./sx).^2+((y-double(muy))./sy).^2+((X1-k*x.*y.^-1)./double(sX)).^2));
drfTHX=integral2(nrfTHXh,xmin,xmax,ymin,ymax);
d=1/drfTHX;
fTHX=d*nrfTHX;
mDx=vpa(int(nrfTHX,y,ymin,ymax),5);
xhati=int(x*mDx,x,xmin,xmax);
xhat=d*vpa(xhati,5);
g0=nrfTHXh(mux,muy);g1x=nrfTHXh(x,muy);g1y=nrfTHXh(mux,y);
mDx1=vpa(int(g1x+g1y-g0,y,ymin,ymax),5);
xhat1=vpa(int(x*mDx1,x,xmin,xmax),5);
xMSEi=int((x-xhat)^2*mDx,x,xmin,xmax);
xMSE0=d*vpa(xMSEi,5);
% mDy=int(nrfTHX,x,xmin,xmax);
% yhati=int(y*mDy,y,ymin,ymax);
% yhat=vpa(yhati,5);
mDy=vpa(int(g1x+g1y-g0,x,xmin,xmax),5);
yhat=d*vpa(int(y*mDy,y,ymin,ymax),5);
yMSEi=int((y-yhat)^2*mDy,y,ymin,ymax);
yMSE0=d*vpa(yMSEi,5);
xhat0=xhat;yhat0=yhat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muXh=@(x1,x2) k*x1*vpa((int((1/x)*normpdf(x,x2,sigth2),x,th2min,th2max)),5);
for j=1:15
mu1in=xhat;mu2in=yhat;
syms x y
muXin=muXh(mu1in,mu2in);
c1=double(vpa(muXin/k,5));
Ig=@(x,y) ((x/y)-c1).^2*exp(-0.5*((x-mu1in).^2+((y-mu2in)./sigth2).^2));
sig2X=k^2*(1/(2*pi*sigth1*sigth2))*vpa(int(int(Ig(x,y),x,th1min,th1max),y,th2min,th2max),5);
sX=sqrt(sig2X);
syms x y
nrfTHX=exp(-0.5*(((x-mu1in)/sx)^2+((y-mu2in)/sy)^2+((X1-k*x*y^-1)/sX)^2));
nrfTHXh=@(x,y) exp(-0.5.*(((x-double(mu1in))./sx).^2+((y-double(mu2in))./sy).^2+((X1-k*x.*y.^-1)./double(sX)).^2));
drfTHX=integral2(nrfTHXh,xmin,xmax,ymin,ymax);
d=1/drfTHX;
fTHX=d*nrfTHX;
mDx=vpa(int(nrfTHX,y,ymin,ymax),5);
xhati=int(x*mDx,x,xmin,xmax);
xhat=d*vpa(xhati,5);
xhatM(j)=xhat;
g0=nrfTHXh(mu1in,mu2in);g1x=nrfTHXh(x,mu2in);g1y=nrfTHXh(mu1in,y);
mDx1=vpa(int(g1x+g1y-g0,y,ymin,ymax),5);
xhat1=vpa(int(x*mDx1,x,xmin,xmax),5);
xMSEi=int((x-xhat)^2*mDx,x,xmin,xmax);
xMSE=d*vpa(xMSEi,5);
xMSEM(j)=xMSE;
mDy=vpa(int(g1x+g1y-g0,x,xmin,xmax),5);
yhat=d*vpa(int(y*mDy,y,ymin,ymax),5);
yhatM(j)=yhat;
yMSEi=int((y-yhat)^2*mDy,y,ymin,ymax);
yMSE=d*vpa(yMSEi,5);
yMSEM(j)=yMSE;
end
figure(1)
scatter(1:16,[xhat0,xhatM])
set(gca,'fontsize',10)
title ('Estimated values of load P by Bayesian Inference with HDMR','FontSize',12,'FontWeight','bold')
xlabel ('No. of iterations','FontSize',12,'FontWeight','bold')
ylabel ('Estimated P','FontSize',12,'FontWeight','bold')
grid on
figure(2)
scatter(1:16,[xMSE0,xMSEM])
set(gca,'fontsize',10)
title ('Mean Squared Error for P by Bayesian Inference with HDMR','FontSize',12,'FontWeight','bold')
xlabel ('No. of iterations','FontSize',12,'FontWeight','bold')
ylabel ('MSE values','FontSize',12,'FontWeight','bold')
grid on
figure (3)
scatter(1:16,[yhat0,yhatM])
set(gca,'fontsize',10)
title ('Estimated values of flexural rigidity EI by Bayesian Inference with HDMR','FontSize',12,'FontWeight','bold')
xlabel ('No. of iterations','FontSize',12,'FontWeight','bold')
ylabel ('Estimated EI values','FontSize',12,'FontWeight','bold')
grid on
figure (4)
scatter(1:16,[yMSE0,yMSEM])
set(gca,'fontsize',10)
title ('Mean Squared Error for flexural rigidity EI by Bayesian Inference with HDMR','FontSize',12,'FontWeight','bold')
xlabel ('No. of iterations','FontSize',12,'FontWeight','bold')
ylabel ('MSE values','FontSize',12,'FontWeight','bold')
grid on 