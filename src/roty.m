%ROTY Rotation about Y axis
%
% R = ROTY(THETA) is a rotation matrix representing a rotation of THETA 
% radians about the y-axis.
%
% R = ROTY(THETA, 'deg') as above but THETA is in degrees.
%
% See also ROTX, ROTZ, ANGVEC2R.


% Copyright (C) 1993-2011, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for Matlab (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.

function R = roty(t, deg)
    if nargin > 1 && strcmp(deg, 'deg')
        t = t *pi/180;
    end
    ct = cos(t);
    st = sin(t);
    R = [
        ct  0   st 0
        0   1   0 0 
       -st  0   ct 0
       0 0 0 1
       ];
