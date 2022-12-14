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

M = M_test(ni,dx,dt,au,ap,f);
Madj = M';
% get trajectory from initial state iterate
x_prop(:,1) = x0;
for k = 1:1:kend-1
    x_prop(:,k+1) = M*x_prop(:,k);
end

%  add obs contribution to gradient
for k = kend:-1:1
   if (mod(k,obs_freq)==0) % if obs exists at timestep k
      innov = xobs(:,k) - H*x_prop(:,k);
      grd = grd + H'*Rinv*innov;
   end
%  --- propagate gradient of obs term at time t backward one timestep ---
   grd = Madj*grd;
end

% add contribution from background
grd = -grd + Binv*bdif;

