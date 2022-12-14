function Tn = upwindad (Tnp1, t)

%UPWINDAD  ADJOINT of UPWIND difference scheme 
%
%   Global parameters:
%       Courant_Number:  u  * dt / dx
%       Diffusivity:     mu * dt / (dx^2)
%
%  where u is the advective speed; mu is the physical diffusivity
%  parameter; dx is the grip point distance; and dt is the time
%  interval.
%
%  Notice that stability requires:  
%
%      dt <=  0.5 * dx^2 / ( 1 + 0.5 * R )
%
%  where R = Courant_Number / Diffusivity
%
%Saroja Polavarapu, 20 Feb 01  Initial code.

global  Courant_Number Diffusivity

C = Courant_Number;
s = Diffusivity;
jdim = size(Tnp1,1);
Tn = zeros(jdim,1);

% --- Start reversing code ---

%
%  Adjoint of Update periodic boundaries
%
  j = jdim;
% forward: Tnp1(j)  = (s+C)*Tn(j-1)  + (1-2*s-C)*Tn(j) + s*Tn(1); 
  Tn(1)   = Tn(1)   +         s*Tnp1(j);
  Tn(j)   = Tn(j)   + (1-2*s-C)*Tnp1(j);
  Tn(j-1) = Tn(j-1) +     (s+C)*Tnp1(j);
% Tnp1(j) = 0.0;

  j = 1;
% forward: Tnp1(j)  = (s+C)*Tn(jdim) + (1-2*s-C)*Tn(j) + s*Tn(j+1);
  Tn(j+1)  = Tn(j+1)  +         s*Tnp1(j);
  Tn(j)    = Tn(j)    + (1-2*s-C)*Tnp1(j);
  Tn(jdim) = Tn(jdim) +     (s+C)*Tnp1(j);
% Tnp1(j)  = 0.0;
%
%  Adjoint of Loop over the central domain
%
for j = jdim-1:-1:2
% forward: Tnp1(j) = C*Tn(j-1) + (1-C)*Tn(j);
  Tn(j)   = Tn(j)   + (1-C)*Tnp1(j);
  Tn(j-1) = Tn(j-1) +     C*Tnp1(j);
% Tnp1(j) = 0.0;
end 
