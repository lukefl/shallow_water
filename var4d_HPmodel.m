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
dt=1.0*3600/nk;        %Time step in s

% Assimilation parameters
obs_freq = 1;
obs_hole = false;
std_obs = 0.02; % relative error std dev
std_bkg = 1.0; % relative error std dev
Ld = ni*dx/4; % background error correlation scales

% matrix to scale relative errors for covariance matrices
scalingmat = zeros(3*ni,3*ni);
for i=1:2*ni
    for j=1:2*ni
        scalingmat(i,j) = au;
    end
    for j=2*ni:3*ni
        scalingmat(i,j) = sqrt(au*ap);
    end
end
for i=2*ni:3*ni
    for j=1:2*ni
        scalingmat(i,j) = sqrt(au*ap);
    end
    for j=2*ni:3*ni
        scalingmat(i,j) = ap;
    end
end
scalingmat_sq = scalingmat.*scalingmat;

% initialize assimilation variables
% background error
Cor0 = gcorr('gauss', ni*dx, Ld/100, 3*ni, 0);
Cor0 = Cor0/Cor0(1,1);
% B = (std_bkg^2)*scalingmat_sq.*Cor0;
% Bsqrt = std_bkg*scalingmat.*Cor0;
B = (std_bkg^2)*scalingmat_sq.*eye(3*ni);
Bsqrt = std_bkg*scalingmat.*eye(3*ni);
Binv = inv(B);
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
% observation error
R = std_obs^2*scalingmat_sq*eye(3*nobs);
Rsqrt = std_obs*scalingmat.*eye(3*nobs);
R = eye(3*nobs);
Rinv = inv(R);

% figure
% ax = axes;
% contourf(log10(abs(B)))
% colorbar()
% set(ax,'Ydir', 'reverse');
% % set(ax,'Xdir', 'reverse');

%Initialize random number generator to ensure the same
%sequence of random numbers each time the model is run
rand('state',0);

%Set initial conditions
[u1,v1,p1]=set_init(ni,ap,dx,ra,f,1);

% Define initial state vector
xold=(set_state(u1,v1,p1,ni))';

x(:,1)=xold;
xobs=zeros(3*nobs,nk);

M = M_test(ni,dx,dt,au,ap,f);
M_ad = M'; % adjoint
% x_test = zeros(3*ni,nk);
% x_test(:,1) = xold;
%%%%%%%%%%%%%%%% loop in time %%%%%%%%%%%%%
for k = 1:kend-1
	xold = x(:,k);
	xobs(:,k) = H*xold + 0*Rsqrt*randn(3*nobs,1);
	xnew = (HP_solver(xold', ni, dx, au, ap, f, ra, k, dt, nfor))';
% 	x(:,k+1) = xnew;
	% test M matrix
    x(:,k+1) = M*x(:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(num2str((x_test(kend)'*x_test(kend)-x(kend)'*x(kend))/(x(kend)'*x(kend)),15))
% 
% [uf,vf,pf]=get_uvp((x(:,kend))',ni);
% [uf_test,vf_test,pf_test]=get_uvp((x_test(:,kend))',ni);
% [u1,v1,p1]=get_uvp((x(:,1))',ni);
% [u1_test,v1_test,p1_test]=get_uvp((x_test(:,kend))',ni);
% 
% % Plot initial condition
% xgrid=(1:ni)*dx*0.001;
% xgrid_obs=(1:nobs)*dx*0.001;
% figure(1),subplot(3,1,1);
% plot(xgrid,vf,'k-',xgrid,vf_test,'r-'),title('Final conditions'),
% ylabel('u (m/s)'),legend('old method','matrix method');
% subplot(3,1,2);plot(xgrid,u1,'k-'),ylabel('v (m/s)');
% subplot(3,1,3);plot(xgrid,v1,'k-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

xback=zeros(3*ni,nk);
xback(:,kend)=x(:,kend);
%%%%%%%%%%%%%%%% loop in backward in time %%%%%%%%%%%%%
for k = kend:-1:2
	xold = xback(:,k);
% 	xnew = (HP_adjoint(xold', ni, dx, au, ap, f, ra, k, dt, nfor))';
    xnew = M_ad*xold;
	xback(:,k-1) = xnew;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u1b,v1b,p1b]=get_uvp((xback(:,1))',ni);


%%%%%%%%%%%%%%% 4D-Var section %%%%%%%%%%%%%
% Cost function
xb0 = zeros(3*ni,1); % background for time=0
% xb0 = x(:,1)+10*ones(3*ni,1);
% xobs = x;
% xobs = zeros(3*nobs,nk);
% xtest = xb0;
% xtest = xobs(:,1);
xtest = ones(3*ni,1);

Jcost = cost(xtest,xobs,xb0,Rinv,Binv,H,dt);
grd = gradcost(xtest,xobs,xb0,Rinv,Binv,H,dt);
grdnorm = grd'*grd

%   Approximation to the Hessian matrix
psi = eye(3*ni);
alpha = 0.1^3;
for i = 1:3*ni
	pert = psi(:,i);
	xp = xtest + alpha*pert;
	costp = cost(xp,xobs,xb0,Rinv,Binv,H,dt);
	grdp = gradcost(xp,xobs,xb0,Rinv,Binv,H,dt);
	Hess(:,i) = (grdp - grd)/alpha;
end

%   Gradient test
Hesschk = 0.5*grd'*Hess*grd;
for loop = 1:1:10
	alpha = 0.1^(loop+2);
	xp = xtest + alpha*grd;
	costp = cost(xp,xobs,xb0,Rinv,Binv,H,dt);
	Jdifp = (costp-Jcost)/(alpha*grdnorm);
	ratio2 = (costp-Jcost-alpha*grdnorm)/(alpha^2*Hesschk);
disp (num2str([alpha Jdifp ratio2],15));
end

%   Solution using Newton's method
xa = zeros(3*ni,nk);
Hessi = inv(Hess);
xa(:,1) = xb0 - Hessi*grd; % xa = analysis
xa0 = xa(:,1); % analysis at t=0

%   Integrate solution (at initial time) to get final time
for k=1:kend-1
%     xa(:,k) = HP_solver(xa(:,k-1), ni, dx, au, ap, f, ra, k, dt, nfor);
    xa(:,k+1) = M*xa(:,k);
end
[u1a,v1a,p1a]=get_uvp((xback(:,kend))',ni);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot initial condition
xgrid=(1:ni)*dx*0.001;
xgrid_obs=(1:nobs)*dx*0.001;
figure(1),subplot(3,1,1);
plot(xgrid,u1,'k-',xgrid,u1a,'r-'),title('Initial conditions'),
ylabel('u (m/s)'),legend('true state','analysis');
subplot(3,1,2);plot(xgrid,v1,'k-'),ylabel('v (m/s)');
subplot(3,1,3);plot(xgrid,p1,'k-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

% Define final conditions
[u1,v1,p1]=get_uvp(x(:,kend)',ni);
[u1_obs,v1_obs,p1_obs]=get_uvp(xobs(:,kend)',nobs);

% Plot final conditions
figure(2),subplot(3,1,1);
plot(xgrid,u1,'k-',xgrid,u1a),title('Final conditions'),
ylabel('u (m/s)'),legend('reference state','perturbed state');
subplot(3,1,2);plot(xgrid,v1,'k-'),ylabel('v (m/s)');
subplot(3,1,3);plot(xgrid,p1,'k-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

% Plot time series at midpoint of domain
t=(1:kend+1)*dt/3600.0;
ut1=x(imid,:);
vt1=x(ni+imid,:);
pt1=x(2*ni+imid,:);

figure(3),subplot(3,1,1);
plot(t,ut1,'k-'),title('Time series'),
ylabel('u (m/s)'),legend('reference state','perturbed state');
subplot(3,1,2);plot(t,vt1,'k-'),ylabel('v (m/s)');
subplot(3,1,3);plot(t,ut1.^2+vt1.^2),xlabel('t (hr)'),ylabel('p (m^2/s^2)');

% while(waitforbuttonpress()==0) pause(1) end