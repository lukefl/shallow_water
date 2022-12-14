%
% One-Dimensional linearized shallow water model in a periodic domain
%
% Hinkelmann-Phillips model in 1D
% Reference: Daley (1991) ch. 6.4
% Saroja Polavarapu Jan. 14, 2004

%Initialize workspace
clear all;
clf,figure(1);
clf,figure(2);
clf,figure(3);

%first set the model grid parameters.
nk=360;                 %Number of time steps
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
dt=72.0*3600/nk;        %Time step in s

%Initialize random number generator to ensure the same
%sequence of random numbers each time the model is run
rand('state',0);

%Set initial conditions
[u1,v1,p1]=set_init(ni,ap,dx,ra,f,3);
[u2,v2,p2]=set_init(ni,ap,dx,ra,f,3);
 
% Plot initial condition
xgrid=(1:ni)*dx*0.001;
figure(1),subplot(3,1,1);
plot(xgrid,u1,'k-',xgrid,u2,'r-'),title('Initial conditions'),
ylabel('u (m/s)'),legend('reference state','perturbed state');
subplot(3,1,2);plot(xgrid,v1,'k-',xgrid,v2,'r-'),ylabel('v (m/s)');
subplot(3,1,3);plot(xgrid,p1,'k-',xgrid,p2,'r-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

% Define initial state vector
xold=set_state(u1,v1,p1,ni);
yold=set_state(u2,v2,p2,ni);
      
x(1,:)=xold;
y(1,:)=yold;

%%%%%%%%%%%%%%%% loop in time %%%%%%%%%%%%%
for k = 1:kend;
	xold = x(k,:);
   yold = y(k,:);
   xnew = HP_solver(xold, ni, dx, au, ap, f, ra, k, dt, nfor);
	ynew = HP_solver(yold, ni, dx, au, ap, f, ra, k, dt, nfor);
	x(k+1,:) = xnew;
	y(k+1,:) = ynew;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define final conditions
[u1,v1,p1]=get_uvp(xnew,ni);
[u2,v2,p2]=get_uvp(ynew,ni);

% Plot final conditions
figure(2),subplot(3,1,1);
plot(xgrid,u1,'k-',xgrid,u2,'r-'),title('Final conditions'),
ylabel('u (m/s)'),legend('reference state','perturbed state');
subplot(3,1,2);plot(xgrid,v1,'k-',xgrid,v2,'r-'),ylabel('v (m/s)');
subplot(3,1,3);plot(xgrid,p1,'k-',xgrid,p2,'r-'),xlabel('x (km)'),ylabel('p (m^2/s^2)');

% Plot time series at midpoint of domain
t=(1:kend+1)*dt/3600.0;
ut1=x(:,imid);
vt1=x(:,ni+imid);
pt1=x(:,2*ni+imid);
ut2=y(:,imid);
vt2=y(:,ni+imid);
pt2=y(:,2*ni+imid);

figure(3),subplot(3,1,1);
plot(t,ut1,'k-',t,ut2,'r-'),title('Time series'),
ylabel('u (m/s)'),legend('reference state','perturbed state');
subplot(3,1,2);plot(t,vt1,'k-',t,vt2,'r-'),ylabel('v (m/s)');
subplot(3,1,3);plot(t,ut1.^2+vt1.^2,'k-',t,pt2,'r-'),xlabel('t (hr)'),ylabel('p (m^2/s^2)');

