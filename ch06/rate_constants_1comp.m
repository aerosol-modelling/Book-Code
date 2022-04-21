function [K,E,S] = rate_constants_1comp(nmax)

% This function calculates the rate constants for simulating the cluster concentrations
% Output:
%    K = collision rate coefficient matrix (m^3 s^-1)
%    E = evaporation rate coefficient matrix (s^-1)
%    S = cluster scavenging coefficient vector (s^-1)

% Physical constants
kB = 1.3806504e-23;          % Boltzmann constant (J K^-1)
amu2kg = 1e-3/6.02214179e23; % Conversion from atomic mass units to kg

% Temperature
temp = 280;                  % (K)

% Cluster properties
m1 = 100*amu2kg;             % Mass of one molecule (kg)
rho = 1500;                  % Assumed density (kg m^-3)

% Cluster mass, volume and diameter
m = m1.*(1:nmax);            % (kg)
V = m./rho;                  % (m^3)
dp = (6/pi.*V).^(1/3);       % (m)

% Allocate the matrices and vectors for the rate constants
K = zeros(nmax,nmax);
E = zeros(nmax,nmax);
S = zeros(1,nmax);

for i = 1:nmax
    for j = i:nmax
        % Hard-sphere collision rate of clusters i and j
        K(i,j) = (3/4/pi)^(1/6)*(6*kB*temp*(1/m(i)+1/m(j)))^(1/2)*(V(i)^(1/3)+V(j)^(1/3))^2;
        K(j,i) = K(i,j);
    end
    % Loss rate constant to external sinks
    S(i) = 1e-3*(dp(i)/dp(1))^-1.6;
end

% Evaporation rate of molecules from clusters (no fissions assumed)
% A hypothetical example where the smallest clusters evaporate at the order of magnitude of
% quantum-chemistry-based rates for H2SO4-NH3 clusters
% NOTE: When calculated properly, E is dependent on T and also on K due to the detailed balance
E(1,1) = 2e-3;
E(2,1) = 6e-3;
E(3,1) = 3e-3;
E(4,1) = 2e-4;
% Assume a low evaporation rate for the larger sizes
E(5:nmax,1) = 1e-5;
for i = 1:nmax
    E(1,i) = E(i,1);
end

% Divide the constants by two for collisions and evaporations of two identical clusters to
% avoid counting the same process twice
for i = 1:nmax
    K(i,i) = 0.5*K(i,i);
    E(i,i) = 0.5*E(i,i);
end

end