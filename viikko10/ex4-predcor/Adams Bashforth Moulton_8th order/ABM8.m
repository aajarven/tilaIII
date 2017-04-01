%--------------------------------------------------------------------------
% 
% Adams-Bashforth-Moulton 8th-order
% 
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function Y = ABM8(func, f_hist, h)

bash = [434241.0, -1152169.0, 2183877.0, -2664477.0, 2102243.0, -1041723.0, 295767.0, -36799.0];
moul = [36799.0, 139849.0, -121797.0, 123133.0, -88547.0, 41499.0, -11351.0, 1375.0];
divisor = 1.0/120960.0;

% Calculate the predictor using the Adams-Bashforth formula 
Y = f_hist(8,2:7)' + h*divisor* ...
    ( bash(1)*func(f_hist(8,1),f_hist(8,2:7)') + bash(2)*func(f_hist(7,1),f_hist(7,2:7)') + ...
      bash(3)*func(f_hist(6,1),f_hist(6,2:7)') + bash(4)*func(f_hist(5,1),f_hist(5,2:7)') + ...
      bash(5)*func(f_hist(4,1),f_hist(4,2:7)') + bash(6)*func(f_hist(3,1),f_hist(3,2:7)') + ...
      bash(7)*func(f_hist(2,1),f_hist(2,2:7)') + bash(8)*func(f_hist(1,1),f_hist(1,2:7)') );

% Calculate the corrector using the Adams-Moulton formula
Y = f_hist(8,2:7)' + h*divisor* ...
    ( moul(1)*func(f_hist(8,1)+h,Y) + moul(2)*func(f_hist(8,1),f_hist(8,2:7)') + ...
      moul(3)*func(f_hist(7,1),f_hist(7,2:7)') + moul(4)*func(f_hist(6,1),f_hist(6,2:7)') + ...
      moul(5)*func(f_hist(5,1),f_hist(5,2:7)') + moul(6)*func(f_hist(4,1),f_hist(4,2:7)') + ...
      moul(7)*func(f_hist(3,1),f_hist(3,2:7)') + moul(8)*func(f_hist(2,1),f_hist(2,2:7)') );

