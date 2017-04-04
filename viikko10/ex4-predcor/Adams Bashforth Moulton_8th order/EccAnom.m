%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the eccentric anomaly for elliptic orbits
%
% Inputs:
%   M         Mean anomaly in [rad]
%   e         Eccentricity of the orbit [0,1]
% 
% Outputs:
%             Eccentric anomaly in [rad]
%
%--------------------------------------------------------------------------
function [E]  = EccAnom (M, e)

maxit = 15;
i = 1;

% Starting value
M = mod(M, 2.0*pi);

if (e<0.8)
    E=M; 
else
    E=pi;
end
  
  f = E - e*sin(E) - M;
  E = E - f / ( 1.0 - e*cos(E) );
  
% Iteration  
while (abs(f) > eps)

  f = E - e*sin(E) - M;
  E = E - f / ( 1.0 - e*cos(E) );
  i = i+1;
  if (i==maxit)
     error(' convergence problems in EccAnom');
  end
  
end

