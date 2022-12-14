%
% One-Dimensional linearized shallow water model in a periodic domain
%
% Hinkelmann-Phillips model in 1D
% Reference: Daley (1991) ch. 6.4
% Saroja Polavarapu Jan. 14, 2004

%Initialize workspace
echo off
clear all;
%clf,figure(1);
%clf,figure(2);
%clf,figure(3);

%first set the model grid parameters.
% set physical constants and numerical parameters to global
global ra f ap au
global nk ni dx dt obs_freq nobs kend nfor
nk=60;                 %Number of time steps
ni=80;                 %No. of grid points in the zonal direction
imid=(ni)/2+1;         %midpoint of the domain
ni3=ni*3.0;				  %Size of state vector, X
nfor=0;                %Switch to turn on forcing (nfor=0 or 1)
kend=nk;               %Length of integration in timesteps

%Physical constants
ra=6.4e6;              %Radius of the Earth in m
f=1.0e-4;              %Coriolis parameter in 1/s
ap=1.0e5;              %phi_0 scaling factor m^2/s^2
au=20.0;               %U = mean zonal wind in m/s

%Grid calculations
dx=2.0*pi*ra/ni;       %grid spacing in m
dt=6.0*3600/nk;        %Time step in s

% Assimilation parameters
obs_freq = 1;
obs_hole = false;
std_obs = 0.02; % relative error std dev
std_bkg = 1.0; % relative error std dev
Ld = ni*dx/40; % background error correlation scale

% observation network
if (obs_hole==false)
	nobs = ni;
	H = eye(3*nobs);
elseif (obs_hole==true)
	nobs = ni/2;
	H = zeros(3*nobs,3*ni);
	for i=1:1:3*nobs
		H(i,i) = 1.0;
    end
end

% matrix to scale relative errors for covariance matrices
cor1 = gcorr('gauss', ni*dx, Ld, ni, 0);
cor1 = cor1/cor1(1,1);
cor1 = [cor1 zeros(ni, 2*ni)];
cor1 = [cor1; zeros(2*ni, 3*ni)];
cor2 = circshift(cor1,[ni ni]);
cor3 = circshift(cor1,[2*ni 2*ni]);
B = au^2*std_bkg^2*cor1 + au^2*std_bkg^2*cor2 + ap^2*std_bkg^2*cor3;
Binv = inv(B);
Bsqrt = B^(1/2);
% figure(0)
% imshow(B, [0 1])
% colorbar()

% initialize assimilation variables
% observation error
I1 = eye(nobs);
I1 = [I1 zeros(nobs, 2*nobs)];
I1 = [I1; zeros(2*nobs, 3*nobs)];
I2 = circshift(I1, [nobs nobs]);
I3 = circshift(I1, [2*nobs 2*nobs]);
R = au^2*std_obs^2*I1 + au^2*std_obs^2*I2 + ap^2*std_obs^2*I3;
Rsqrt = R^(1/2);
Rinv = inv(R);
% figure(10)
% imshow(log10(R), [0 5])
% colorbar()

%Initialize random number generator to ensure the same
%sequence of random numbers each time the model is run
rand('state',0);

%Set initial conditions
[u0,v0,p0]=set_init(ni,ap,dx,ra,f,1);

% Define initial state vector
xold=(set_state(u0,v0,p0,ni))';

x(:,1)=xold;
xobs=zeros(3*nobs,nk);

M = M_test(ni,dx,dt,au,ap,f);
M_ad = M'; % adjoint
%%%%%%%%%%%%%%%% loop in time %%%%%%%%%%%%%
for k = 1:kend-1
    if mod(k,obs_freq)==0
	    xobs(:,k) = H*x(:,k) + 0*Rsqrt*randn(3*nobs,1);
    end
    x(:,k+1) = M*x(:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% 4D-Var section %%%%%%%%%%%%%
% Cost function
xb0 = x(:,1) + Bsqrt*randn(3*ni,1); % background for time=0
x_guess = xb0; %zeros(3*ni,1);

Jcost = cost(x_guess,xobs,xb0,Rinv,Binv,H,dt);
grd = gradcost(x_guess,xobs,xb0,Rinv,Binv,H,dt);
grdnorm = grd'*grd;

%   Approximation to the Hessian matrix
psi = eye(3*ni);
alpha = 0.1^3;
for i = 1:3*ni
	pert = psi(:,i);
	xp = x_guess + alpha*pert;
	costp = cost(xp,xobs,xb0,Rinv,Binv,H,dt);
	grdp = gradcost(xp,xobs,xb0,Rinv,Binv,H,dt);
	Hess(:,i) = (grdp - grd)/alpha;
end

%   Gradient test
Hesschk = 0.5*grd'*Hess*grd;
for loop = 1:1:5
	alpha = 0.1^(loop+2);
	xp = x_guess + alpha*grd;
	costp = cost(xp,xobs,xb0,Rinv,Binv,H,dt);
	Jdifp = (costp-Jcost)/(alpha*grdnorm);
	ratio2 = (costp-Jcost-alpha*grdnorm)/(alpha^2*Hesschk);
disp (num2str([alpha Jdifp ratio2],15));
end

%   Solution using Newton's method
xa = zeros(3*ni,nk);
Hessi = inv(Hess);
xa(:,1) = x_guess - Hessi*grd; % xa = analysis
xa0 = xa(:,1); % analysis at t=0

%   Integrate solution (at initial time) to get final time
for k=1:kend-1
    xa(:,k+1) = M*xa(:,k);
end

% Get fields to plot
[u0b,v0b,p0b]=get_uvp(xb0',ni);
[uf,vf,pf]=get_uvp(x(:,kend)',ni);
[uf_obs,vf_obs,pf_obs]=get_uvp(xobs(:,kend)',nobs);
[ufa,vfa,pfa]=get_uvp((xa(:,kend))',ni);
[u0a,v0a,p0a]=get_uvp((xa(:,1))',ni);

% Plot initial condition
xgrid=(1:ni)*dx*0.001;
xgrid_obs=(1:nobs)*dx*0.001;
figure(1),subplot(3,1,1);
plot(xgrid,u0,'k-',xgrid,u0b,'b-',xgrid,u0a,'r-'),title('Initial conditions'),
ylabel('u (m/s)'),legend('true state','background','analysis');
subplot(3,1,2);plot(xgrid,v0,'k-',xgrid,v0b,'b-',xgrid,v0a,'r-'),ylabel('v (m/s)');
subplot(3,1,3);plot(xgrid,p0,'k-',xgrid,p0b,'b-',xgrid,p0a,'r-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

% Plot final conditions
figure(2),subplot(3,1,1);
plot(xgrid,uf,'k-',xgrid,ufa,'r-'),title('Final conditions'),
ylabel('u (m/s)'),legend('true state','analysis');
subplot(3,1,2);plot(xgrid,vf,'k-',xgrid,vfa,'r-'),ylabel('v (m/s)');
subplot(3,1,3);plot(xgrid,pf,'k-',xgrid,pfa,'r-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

% Plot time series at midpoint of domain
t=(1:kend+1)*dt/3600.0;
ut1=x(imid,:);
vt1=x(ni+imid,:);
pt1=x(2*ni+imid,:);

% figure(3),subplot(3,1,1);
% plot(t,ut1,'k-'),title('Time series'),
% ylabel('u (m/s)'),legend('reference state','perturbed state');
% subplot(3,1,2);plot(t,vt1,'k-'),ylabel('v (m/s)');
% subplot(3,1,3);plot(t,ut1.^2+vt1.^2),xlabel('t (hr)'),ylabel('p (m^2/s^2)');

% while(waitforbuttonpress()==0) pause(1) end