function x=set_state(u,v,p,ni)
%  Define state vector, x = [u,v,p]
x(1:ni)=u;
x(ni+1:2*ni)=v;
x(2*ni+1:3*ni)=p;
