%       TESTADJ
%       tests adjoint code for HP_solver.m (HP_adjoint.m)

%first set the model grid parameters.
nk=30;                 %Number of time steps
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
dt=6*3600/nk;        %Time step in s

% Assimilation parameters
obs_freq = 1;
obs_hole = true;
std_obs = 5;

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
%Initialize random number generator to ensure the same
%sequence of random numbers each time the model is run
rand('state',0);

%Set initial conditions
[u1,v1,p1]=set_init(ni,ap,dx,ra,f,3);

% Define initial state vector
xi=set_state(u1,v1,p1,ni)';
x=xi;
xim = xi;
xm = xim;
M = M_test(ni,dx,dt,au,ap,f);
%%%%%%%%%%%%%%%% loop in time %%%%%%%%%%%%%
for k = 1:kend
    disp(k)
	x = (HP_solver(x', ni, dx, au, ap, f, ra, k, dt, nfor))';
	xm = M*xm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xf = x;
xmf = xm;

%%%%%%%%%%%%%%%% loop in backward in time %%%%%%%%%%%%%
for k = kend:-1:1
    disp(k)
	xm = M'*xm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xa=xm;

echo off
% --- adjoint test, Problem 3a ---
% for Y = L X the adjoint is X* = L* Y
% The adjoint test checks that <Y,LX> = <L*Y,X>
% where lhs = <L*Y,X> and rhs = <Y,LX> 
lhs = xim'*xa
rhs = xmf'*xmf
disp (num2str([lhs rhs],15));

