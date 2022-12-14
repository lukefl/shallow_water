%--------------------------------------------------------------
function dXdt = HP_rhs(X,ni,dx,au,ap,f,ra,k,dt,nfor)

%Right Hand Side of the Hinkelmann-Phillips model
%Jan. 8 Saroja Polavarapu

amf1=1.0;             %amplitude of large scale forcing
amf2=0.5;             %amplitude of small scale forcing
tforce=6.0*3600.0;    %timescale of forcing in sec

%   Initializae local variables
dXdt = zeros(size(X),size(X)); 
[u,v,p]=get_uvp(X,ni);

%   Take x-derivative of u, v and p using period boundary conditions
dxi2=1.0/(2.0*dx);
for i = 2:ni-1
   dudx(i)=dxi2*(X(i+1)-X(i-1));
   dvdx(i)=dxi2*(X(ni+i+1)-X(ni+i-1));
   dpdx(i)=dxi2*(p(2*ni+i+1)-p(2*ni+i-1));
end
dudx(1)=dxi2*(X(2)-X(ni));
dvdx(1)=dxi2*(X(ni+2)-v(2*ni));
dpdx(1)=dxi2*(p(2*ni+2)-p(3*ni));
dudx(ni)=dxi2*(X(1)-X(ni-1));
dvdx(ni)=dxi2*(X(ni+1)-X(2*ni-1));
dpdx(ni)=dxi2*(p(2*ni+1)-p(3*ni-1));

%   Obtain tendencies
for i = 1:ni
   dudt(i) = -au*dudx(i)+f*v(i)-dpdx(i);
   dvdt(i) = -au*dvdx(i)-f*u(i);
   dpdt(i) = -au*dpdx(i)+f*au*v(i)-ap*dudx(i)+...
             nfor*(0.5*ap*au/(2.0*ra))*(1.0+cos(2*pi*dt*k/tforce))...
                 *(amf1*cos(3.0*dx*i/ra)+amf2*sin(7.0*dx*i/ra));
end

%   Define output variable dXdt
dXdt=set_state(dudt,dvdt,dpdt,ni);
