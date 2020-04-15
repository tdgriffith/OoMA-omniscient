function [y,t] = stringFDTD(T,x,rho,C,td,f_idx)
% [y,t] = STRINGFDTD(T,x,rho,C,td,f)
%   Finite difference, time domain solution to the 1D wave equation that
%   describes transverse motion of a taut string with fixed boundary
%   conditions and (optional) viscous damping. This simulation applies a
%   unit impulse to a location specified using f_idx.
%
%   INPUTS:
%   T       string tension (scalar)
%   x       vector of string spatial coordinates
%   rho     mass per unit length of the string (scalar)
%   C       viscous damping coefficient (scalar)
%   td      simulation duration (seconds, scalar)
%   f_idx   index corresponding to spatial vector x to apply impulsive
%           force
%
%   OUTPUTS:
%   y       output of displacement time histories of the string at the
%           locations in the x vector, where each row is the time history
%           of the corresponding spatial location
%   t       time vector (in seconds) corresponding to the columns of y
%
%   NOTES:
%   (1)     the time step is automatically chosen based on the spatial
%           discretization to ensure numerical stability

v = sqrt(T/rho); % wave speed
dx = x(2)-x(1); % spatial discretization
dt = dx/v; % time discretization to ensure stable numerical solution
t = 0:dt:td; % time vector
nt = length(t);
nx = length(x);

f = zeros(nt,nx); % force
f(3,f_idx) = 1; % impulse applied at specified location on second timestep
y = zeros(nt,nx); % assumes stationary initial conditions

for ct = 2:nt-1 %loop over all time
    for cx = 2:nx-1 %loop over all positions but boundaries
        % without damping (for reference)
        % y(ct+1,cx) = dt^2/rho*(T/dx^2*(y(ct,cx+1)-2*y(ct,cx)+y(ct,cx-1))+f(ct,cx))+2*y(ct,cx)-y(ct-1,cx);
        
        % with damping (use C=0 for no damping)
        y(ct+1,cx) = 1/(rho/dt^2+C/(2*dt))*(T/dx^2*(y(ct,cx+1)-2*y(ct,cx)+y(ct,cx-1))+f(ct,cx)-rho/dt^2*(-2*y(ct,cx)+y(ct-1,cx))+C/(2*dt)*y(ct-1,cx));
    end
    % free-free BCs (for reference)
    % y(ct+1,1) = y(ct+1,2); %left boundary no slope for free-free
    % y(ct+1,nx) = y(ct+1,nx-1); %right boundary no slope for free-free
end

y = y/dt; % scaling due to force applied after first timestep
y = y';


end