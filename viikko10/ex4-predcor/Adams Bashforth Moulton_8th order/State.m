%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the satellite state vector from osculating Keplerian elements 
%   for elliptic orbits
%
% Inputs:
%   gm        Gravitational coefficient
%   Kep       Keplerian elements (a,e,i,Omega,omega,M) with
%               a      Semimajor axis 
%               e      Eccentricity 
%               i      Inclination [rad]
%               Omega  Longitude of the ascending node [rad]
%               omega  Argument of pericenter  [rad]
%               M      Mean anomaly at epoch [rad]
%   dt        Time since epoch
% 
% Output:
%             State vector (x,y,z,vx,vy,vz)
%
%   Reference:
%   O.Montenbruck, E. Gill, "Satellite Orbits - Models, Methods,
%   and Applications", Springer Verlag, Heidelberg, (2005)
%--------------------------------------------------------------------------
function [Y] = State ( gm, Kep, dt )

% Keplerian elements at epoch  
a = Kep(1);  Omega = Kep(4);
e = Kep(2);  omega = Kep(5);
i = Kep(3);  M0    = Kep(6);

% Mean anomaly  
if (dt==0.0)
  M = M0;
else
  n = sqrt (gm/(a*a*a));
  M = M0 +n*dt;
end

% Eccentric anomaly  
E  = EccAnom(M,e);

cosE = cos(E);
sinE = sin(E);

% Perifocal coordinates
fac = sqrt ( (1.0-e)*(1.0+e) );

R = a*(1.0-e*cosE);  % Distance
V = sqrt(gm*a)/R;    % Velocity

r = [ a*(cosE-e), a*fac*sinE , 0.0 ]';
v = [ -V*sinE   , +V*fac*cosE, 0.0 ]';

% Transformation to reference system (Gaussian vectors)
PQW = R_z(-Omega) * R_x(-i) * R_z(-omega);

r = PQW*r;
v = PQW*v;

Y = [r;v];

