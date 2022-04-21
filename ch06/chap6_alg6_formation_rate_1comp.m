% Calculate the formation rate out of the simulation system
J = zeros(length(T),1);     % (m^-3 s^-1)
for i = 1:nmax
    for j = i:nmax
        if i+j > nmax
            J = J + K(j,i)*C(:,i).*C(:,j);
        end
    end
end