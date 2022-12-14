function [u1,v1,p1]=set_init(ni,ap,dx,ra,f,icase)
% Define initial conditions for Hinkelmann-Phillips model
gsd=0.2;
amubd=1.6;
amvbd=1.6;
ampbd=1.6;
ubd=2.*amubd*(rand(1,1)-0.5);
vbd=2.*amvbd*(rand(1,1)-0.5);
pbd=2.*ampbd*(rand(1,1)-0.5);
;
ind=[1:ni];
xx=ind*dx/ra;
if icase == 1
   p1=0.1*ap*( 0.8*sin(2*xx) - 0.2*cos(6*xx) + 0.2*sin(5*xx).*sin(3*xx) );
   v1=0.1*ap*( 1.6*cos(2*xx) + 1.2*sin(6*xx) + 0.8*cos(4*xx).*sin(2*xx) + ...
      0.4*sin(4*xx).*cos(2*xx) + gsd*cos(8*xx).*sin(4*xx) )/(f*ra);
   u1=6.*( sin(8*xx).*cos(3*xx) + 0.2*sin(7*xx) );
elseif icase == 2  
   p1=0.1*ap*(0.8*sin(3*xx+pbd) - 0.2*cos(5*xx) + 0.2*sin(4*xx).*sin(2*xx+pbd));
   v1=0.1*ap*(2.4*cos(3*xx+pbd) + 1.0*sin(6*xx) + 0.8*cos(4*xx).*sin(2*xx+pbd)+...
      0.4*sin(4*xx).*cos(2*xx+pbd) + gsd*cos(4*xx+pbd).*cos(8*xx+pbd))/(f*ra);
   u1=6.0*sin(8*xx+ubd).*cos(5*xx);
elseif icase == 3
   p1=0.1*ap*0.8*cos(2*xx);
   v1=ind*0;
   u1=ind*0;
end
