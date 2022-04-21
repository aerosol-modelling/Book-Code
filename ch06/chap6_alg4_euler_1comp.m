% Set the time for which the cluster dynamics are simulated
% Note: Larger systems might take longer to reach the steady state
t_final = 3600;             % (s)

% Solve the equations by simple Eulerian iteration

% Set the time step for the Euler method - no guarantee if this gives correct results...
t_step = 1;                 % (s)
% Number of time steps
imax = floor(t_final/t_step);

T = zeros(1,imax);
C = zeros(imax,nmax);
C(1,:) = c_init;
for i = 2:imax
    dcdt = dgde_1comp(T(i-1),C(i-1,:),K,E,S,Q);
    C(i,:) = C(i-1,:) + dcdt*t_step;
    T(i) = T(i-1) + t_step;
end