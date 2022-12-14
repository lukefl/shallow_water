function Jcost = cost(x0,xobs,xb0,Rinv,Binv,H,dt)
%
%      Input:
%      tobs - observation frequency
%      nobs - no of observation
%      xobs - observations = truth + random error
%      H - observation operator
%      xb - initial background state
%      Binv - inverse of background error cov matrix
%      rr - obs error variance
%     
%      Output:
%      cost - cost function

global ra f ap au nfor 
global ni dx obs_freq kend
bdif  = x0-xb0;
Jbkg = 0.5 * bdif' * Binv * bdif;
Jcost = Jbkg;
x = x0;
M = M_test(ni,dx,dt,au,ap,f);

for k = 1:kend
%  --- If obs available, add contribution to cost ---
   if ( mod(k,obs_freq) == 0)
      xo = xobs(:,k);
      innov = xo-H*x;
      Jobs = 0.5 * innov' * Rinv * innov;
      Jcost = Jcost + Jobs;
   end  
   % update state
   x = M*x;
end   %  End Main Time Loop

