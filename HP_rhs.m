%--------------------------------------------------------------
function dXdt = HP_rhs(X,ni,dx,au,ap,f,ra,k,dt,nfor)

%Right Hand Side of the Hinkelmann-Phillips model
%Jan. 8 Saroja Polavarapu

amf1=1.0;             %amplitude of large scale forcing
amf2=0.5;             %amplitude of small scale forcing
tforce=6.0*3600.0;    %timescale of forcing in sec

%   Initializae local variables
dXdt = zeros(size(X)); 
[u,v,p]=get_uvp(X,ni);

%   Take x-derivative of u, v and p using period boundary conditions
dxi2=1.0/(2.0*dx);
for i = 2:ni-1
   dudx(i)=dxi2*(u(i+1)-u(i-1));
   dvdx(i)=dxi2*(v(i+1)-v(i-1));
   dpdx(i)=dxi2*(p(i+1)-p(i-1));
end
dudx(1)=dxi2*(u(2)-u(ni));
dvdx(1)=dxi2*(v(2)-v(ni));
dpdx(1)=dxi2*(p(2)-p(ni));
dudx(ni)=dxi2*(u(1)-u(ni-1));
dvdx(ni)=dxi2*(v(1)-v(ni-1));
dpdx(ni)=dxi2*(p(1)-p(ni-1));

%   Obtain tendencies
for i = 1:ni
   dudt(i) = -au*dudx(i)+f*v(i)-dpdx(i);
   dvdt(i) = -au*dvdx(i)-f*u(i);
   dpdt(i) = -au*dpdx(i)+f*au*v(i)-ap*dudx(i);
end

%   Define output variable dXdt
dXdt=set_state(dudt,dvdt,dpdt,ni);

    