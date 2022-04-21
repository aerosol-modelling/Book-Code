clear variables

% Set the number of different cluster sizes (i.e. the number of molecules in the largest cluster)
nmax = 5;

% Determine the rate constants
[K,E,S] = rate_constants_1comp(nmax);

% Set the vapor monomer source and no sources for other clusters
Q = zeros(1,nmax);
Q(1) = 5e9;                 % (molecules m^-3 s^-1)

% Set the initial cluster concentrations to zero
c_init = zeros(1,nmax);