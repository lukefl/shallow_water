function grd = gradcost(x0,xobs,xbi,Rinv,Binv,H,dt)
%      Input:
%      tobs - observation frequency
%      nobs - no of observation
%      wo - observations = truth + random error
%      H - observation operator
%      wa - initial background state
%      Binv - inverse of background error cov matrix
%      rr - obs error variance
%     
%      Output:
%      cost - cost function

global ra f ap au nfor 
global nk ni dx obs_freq kend
bdif  = x0-xbi;
grd   = zeros(3*ni,1);
x_prop = zeros(3*ni,nk);

% test for one timestep simulation 
% innov = xobs(:,1) - H*x0;
% grd = -H'*Rinv*innov + Binv*bdif;

M = M_test(ni,dx,dt,au,ap,f);
Madj = M';
% get trajectory from initial state iterate
x_prop(:,1) = x0;
for k = 1:1:kend-1
%     x_prop(:,k+1) = (HP_solver(x_prop(:,k)', ni, dx, au, ap, f, ra, k, dt, nfor))';
%    	xold = x(:,k);
% 	xnew = (HP_solver(xold', ni, dx, au, ap, f, ra, k, dt, nfor))';
% 	x(:,k+1) = xnew;
% 	xobs(:,k) = H*xold + 0*Rsqrt*randn(3*nobs,1);
    x_prop(:,k+1) = M*x_prop(:,k);
end

%  add obs contribution to gradient
for k = kend:-1:1
   if (mod(k,obs_freq)==0) % if obs exists at timestep k
      innov = xobs(:,k) - H*x_prop(:,k);
%       disp(num2str([xobs(:,k)'*xobs(:,k) x_prop(:,k)'*x_prop(:,k) innov'*innov],5))
      grd = grd + H'*Rinv*innov;
   end
%  --- propagate gradient of obs term at time t backward one timestep ---
%    grd = (HP_adjoint(grd', ni, dx, au, ap, f, ra, k, dt, nfor))';
   grd = Madj*grd;
end
% add contribution from t=0
% innov = xobs(:,1) - H*x0;
% grd = grd + H'*Rinv*innov;
% % test if adjoint is properly propagating the gradient
% grd0 = grd;
% for k = 1:1:kend-1
%    grd = (HP_solver(grd', ni, dx, au, ap, f, ra, k, dt, nfor))';
% end
% for k = kend:-1:2
% %  --- propagate gradient of obs term at time t backward one timestep ---
%    grd = (HP_adjoint(grd', ni, dx, au, ap, f, ra, k, dt, nfor))';
% end
% disp([grd'*grd grd0'*grd0 (grd'*grd-grd0'*grd0)/(grd0'*grd0)])

% add contribution from background
grd = -grd + Binv*bdif;

