function xout = HP_adjoint(xin, ni, dx, au, ap, f, ra, k, dt, nfor)
    %integrate the Hinkelmann-Phillips model forward using a 4th order RK scheme    
    
    nv = size(xin,1);
    xx = zeros(nv,1);
    x1 = zeros(nv,1);
    x2 = zeros(nv,1);
    x3 = zeros(nv,1);
    x4 = zeros(nv,1);
    fp = zeros(nv,1);
    w1 = 1.0/6.0;
    w2 = 1.0/3.0;
    w3 = 1.0/3.0;
    w4 = 1.0/6.0;
    
    xx = xin;
    fp = HP_rhs(xx,ni,dx,au,ap,f,ra,k,dt,nfor);
    x1 = dt*fp;
    xx = xin - 0.5*x1;
    fp = HP_rhs(xx,ni,dx,au,ap,f,ra,k,dt,nfor);
    x2 = dt*fp;
    xx = xin - 0.5*x2;
    fp = HP_rhs(xx,ni,dx,au,ap,f,ra,k,dt,nfor);
    x3 = dt*fp;
    xx = xin - x3;
    fp = HP_rhs(xx,ni,dx,au,ap,f,ra,k,dt,nfor);
    x4 = dt*fp;
    xout = xin - (w1*x1 + w2*x2 + w3*x3 + w4*x4);
    % xout = zeros(nv,1);