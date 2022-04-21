clear variables

% Define vapor concentrations (cm^-3)
% NOTE: H2SO4 vapor concentration is here assumed to equal the monomer concentration - this may not
% hold for strong H2SO4-base dimer formation
Ca_vector=10.^(6:0.1:8);                 % Logarithmically evenly spaced H2SO4 values
Cb = 1e9;

J = zeros(size(Ca_vector));

% Loop over the H2SO4 concentrations
for nCa = 1:length(Ca_vector)
    
    Ca = Ca_vector(nCa);
    
    % Set the vapor concentrations
    fid=fopen('driver_input.txt','w');
    fprintf(fid,'constant 1A %e\n',Ca);  % Assuming H2SO4 monomers
    fprintf(fid,'constant 1N %e\n',Cb);
    fclose(fid);
    
    % Call the driver
    [C,T,ok,clust,~,~,~,J_out]=driver_acdc(1e5,'Sources_in','driver_input.txt','repeat');
    
    % Check if the steady state has been reached
    if ok ~= 1
        fprintf('Did not reach steady state for Ca = %.2e cm^-3, Cb = %.2e cm^-3\n',Ca,Cb)
        % Check for negative concentrations
        if ok == -1
            fprintf('Negative concentrations for Ca = %.2e cm^-3, Cb = %.2e cm^-3\n',Ca,Cb)
        end
        continue
    end
    
    % Save the formation rate
    J(nCa) = J_out(end);

end

figure(1)
loglog(Ca_vector,J)
hold on; set(gcf,'Color','white')
xlabel('[H_2SO_4] (cm^{-3})')
ylabel('{\itJ} (cm^{-3} s^{-1})')
