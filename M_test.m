 function M = M_test(ni,dx,dt,au,ap,f)

% % numerical and physical parameters
% nk=30;                 %Number of time steps
% ni=80;                 %No. of grid points in the zonal direction
% ra=6.4e6;              %Radius of the Earth in m
% f=1.0e-4;              %Coriolis parameter in 1/s
% ap=1.0e5;              %phi_0 scaling factor m^2/s^2
% au=20.0;               %U = mean zonal wind in m/s
% dx=2.0*pi*ra/ni;       %grid spacing in m
% dt=6.0*3600/nk;        %Time step in s

I1 = zeros(3*ni);
I2 = zeros(3*ni);
I3 = zeros(3*ni);
% off diagonal matrices (Dp1_1(i,j) = I1(i,j+1), Dm1_1(i,j) = I1(i,j-1), etc)
Dp1_1 = zeros(3*ni,3*ni);
Dm1_1 = zeros(3*ni,3*ni);
Dp1_2 = zeros(3*ni,3*ni);
Dm1_2 = zeros(3*ni,3*ni);
Dp1_3 = zeros(3*ni,3*ni);
Dm1_3 = zeros(3*ni,3*ni);
for i=1:ni
    I1(i,i) = 1;
    I2(ni+i, ni+i) = 1;
    I3(2*ni+i, 2*ni+i) = 1;
end
% I = I1 + I2 + I3;
% % diagonal matrices shifted 1 to the right and 1 to the left from identity matrices
% % diagonal sub-matrices
% Dp1_1 = circshift(I1,1,2);
% Dm1_1 = circshift(I1,-1,2);
% Dp1_2 = circshift(I2,1,2);
% Dm1_2 = circshift(I2,-1,2);
% Dp1_3 = circshift(I3,1,2);
% Dm1_3 = circshift(I3,-1,2);
% % full diagonal matrices
% Dp1 = Dp1_1 + Dp1_2 + Dp1_3;
% Dm1 = Dm1_1 + Dm1_2 + Dm1_3;
% matrices that transform one variable to another
v2u = circshift(I1,ni,2);
u2v = circshift(I1,ni,1);
u2phi = circshift(I1,2*ni,1);
phi2u = circshift(I1,2*ni,2);
v2phi = circshift(I2,ni,1);
% derivative matrices
du2u = circshift(I1,1,2)-circshift(I1,-1,2);
% because variables are stacked in a single vector, circshift doesn't
% handle boundaries correctly - have to do them manually
du2u(ni,ni+1) = 0;
du2u(ni,1) = 1;
du2u(1,3*ni) = 0;
du2u(1,ni) = -1;
%
dv2v = circshift(I2,1,2)-circshift(I2,-1,2);
dv2v(2*ni,2*ni+1) = 0;
dv2v(2*ni,ni+1) = 1;
dv2v(ni+1,ni) = 0;
dv2v(ni+1,2*ni) = -1;
%
dphi2phi = circshift(I3,1,2)-circshift(I3,-1,2);
dphi2phi(3*ni,1) = 0;
dphi2phi(3*ni,2*ni+1) = 1;
dphi2phi(2*ni+1,2*ni) = 0;
dphi2phi(2*ni+1,3*ni) = -1;
%
dphi2u = circshift(phi2u,1,2)-circshift(phi2u,-1,2);
dphi2u(ni,1) = 0;
dphi2u(ni,2*ni+1) = 1;
dphi2u(1,2*ni) = 0;
dphi2u(1,3*ni) = -1;
%
du2phi = circshift(u2phi,1,2)-circshift(u2phi,-1,2);
du2phi(3*ni,ni+1) = 0;
du2phi(3*ni,1) = 1;
du2phi(2*ni+1,3*ni) = 0;
du2phi(2*ni+1,ni) = -1;

% tangent linear model matrix, finally! x(k+1) = M*x(k)
M = eye(3*ni) ...
    + dt*(f*v2u - au*du2u/(2*dx) - dphi2u/(2*dx) ... % dudt
    - f*u2v - au*dv2v/(2*dx) ... % dvdt
    + f*au*v2phi - au*dphi2phi/(2*dx) -ap*du2phi/(2*dx)); % dphidt
% 
% % 
% [u,v,p]=set_init(ni,ap,dx,ra,f,1);
% 
% % Define initial state vector
% x0 = (set_state(u,v,p,ni))';
% %   Take x-derivative of u, v and p using period boundary conditions
% dxi2=1.0/(2.0*dx);
% for i = 2:ni-1
%    dudx(i)=dxi2*(u(i+1)-u(i-1));
%    dvdx(i)=dxi2*(v(i+1)-v(i-1));
%    dpdx(i)=dxi2*(p(i+1)-p(i-1));
% end
% dudx(1)=dxi2*(u(2)-u(ni));
% dvdx(1)=dxi2*(v(2)-v(ni));
% dpdx(1)=dxi2*(p(2)-p(ni));
% dudx(ni)=dxi2*(u(1)-u(ni-1));
% dvdx(ni)=dxi2*(v(1)-v(ni-1));
% dpdx(ni)=dxi2*(p(1)-p(ni-1));
% 
% %   Obtain tendencies
% for i = 1:ni
%    dudt(i) = -au*dudx(i)+f*v(i)-dpdx(i);
%    dvdt(i) = -au*dvdx(i)-f*u(i);
%    dpdt(i) = -au*dpdx(i)+f*au*v(i)-ap*dudx(i);
% end
% 
% %   Define output variable dXdt
% dXdt=set_state(dudt,dvdt,dpdt,ni);
% 
% dXdt_t = (M-eye(3*ni))*x0/dt;
% dudt_t = dXdt(1:ni);
% dvdt_t = dXdt(ni+1:2*ni);
% dpdt_t = dXdt(2*ni+1:3*ni);
% 
% xgrid=(1:ni)*dx*0.001;
% Xgrid=(1:3*ni)*dx*0.001;
% figure(1)
% plot(xgrid,dvdt,'k-',xgrid,dvdt_t,'r-'),title('Final conditions'),
% ylabel('u (m/s)'),legend('explicit version','matrix version');
% 
% 
% % imshow(du2phi, [-1,1]);
% % colorbar()
% % x0 = 500;
% % y0 = 100;
% % width = 500;
% % height = 500;
% % set(gcf,'position',[x0,y0,width,height])