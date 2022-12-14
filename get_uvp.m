function [u,v,p]=get_uvp(x,ni)
%  Get u,v,p from state vector, x
u=x(1:ni);
v=x(ni+1:2*ni);
p=x(2*ni+1:3*ni);
