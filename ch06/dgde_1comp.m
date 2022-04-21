function dcdt = dgde_1comp(t,c,K,E,S,Q)

% This function contains the time derivatives dc/dt of cluster concentrations c
% Input:
%    t = time (can be a "dummy" variable)
%    c = cluster number concentration vector (clusters m^-3)
%    K = collision rate coefficient matrix (m^3 s^-1)
%    E = evaporation rate coefficient matrix (s^-1)
%    S = cluster scavenging coefficient vector (s^-1)
%    Q = cluster source rate vector (clusters m^-3 s^-1)

dcdt = zeros(size(c));
nmax = length(c);

for i = 1:nmax
    
    for j = i:nmax
        
        % The collision removes the smaller clusters...
        K_term = K(j,i)*c(i)*c(j);
        dcdt(i) = dcdt(i) - K_term;
        dcdt(j) = dcdt(j) - K_term;
        
        if i+j <= nmax
            
            % ...and creates a larger cluster.
            dcdt(i+j) = dcdt(i+j) + K_term;
            
            % The corresponding evaporation creates the smaller clusters...
            E_term = E(j,i)*c(i+j);
            dcdt(i) = dcdt(i) + E_term;
            dcdt(j) = dcdt(j) + E_term;
            % ...and removes the larger cluster.
            dcdt(i+j) = dcdt(i+j) - E_term;
            
        end
        
    end
    % Add the source
    dcdt(i) = dcdt(i) + Q(i);
    % Subtract the loss
    dcdt(i) = dcdt(i) - S(i)*c(i);
    
end

end
