clear variables

% Plot the ACDC formation rate as a function of vapor concentrations or temperature

% Generate the equations
% Use the --variable_temp option to easily change the temperature after generating the files
commandstr=['perl acdc_2021_03_23.pl',...
' --i cluster_set_file_acdc.inp',...
' --e formation_free_energy_file_acdc.txt',...
' --variable_temp',...
' --use_cs --cs exp_loss --exp_loss_coefficient 0.001 --exp_loss_exponent -1.6'];
lrun=system(commandstr);
if lrun ~= 0
    error('Running the Perl script failed!')
end
rehash pathreset

% Define vapor concentrations (cm^-3) and temperature (K)
% Here, one of these must be a vector (= the x-axis) and the others are scalars
%Ca = 1e6;              % H2SO4 concentration
Ca = 10.^(6:0.1:8);     % (Logarithmically evenly spaced values)
Cb = 1e10;              % NH3 concentration
%Cb = 10.^(8:0.1:10);
temp = 280;             % Temperature
%temp = 260:10:320;
%IPR = 3;               % (Ion production rate in cm^-3 s^-1 for tests where the cluster set includes charged clusters)

% Define which quantity is on the x-axis - note that also some settings below need to be revised accordingly
% For clarity, this is set manually but as an extra exercise you can design an automatic way to select the x-axis quantity :)
x_vector = Ca;

J = zeros(size(x_vector));

% Loop over the x-axis points
for nxval = 1:length(x_vector)
    
    % Set the x-axis quantity
    Ca = x_vector(nxval);
    %Cb = x_vector(nxval);
    %temp = x_vector(nxval);
    
    % Set the vapor concentrations
    fid=fopen('driver_input.txt','w');
    %fprintf(fid,'constant 1A %e\n',Ca);       % (Assuming H2SO4 monomers)
    fprintf(fid,'constant 1A %e -1A1N\n',Ca);  % Assuming also H2SO4 clustered with base
    fprintf(fid,'constant 1N %e\n',Cb);
    %fprintf(fid,'source neg %e\n',IPR);       % (For tests where the cluster set includes charged clusters)
    %fprintf(fid,'source pos %e\n',IPR);
    fclose(fid);
    
    % Call the driver
    [C,T,ok,clust,~,~,~,J_out]=driver_acdc(1e5,'Sources_in','driver_input.txt','repeat','Temperature',temp);

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
    J(nxval) = J_out(end);
    
end

% % Hint: Note that one could also utilize the existing ACDC Matlab routines in https://github.com/tolenius/ACDC, for example:
% for nxval = 1:length(x_vector)
%     temp = x_vector(nxval);
%     % The wanted output can be directly returned by adding the output variables, here J, in file run_steadystate_ABB.m
%     J(nxval) = run_steadystate_ABB('input_run_steadystate_AN.m','temp',temp);
% end


figure(1)

set(gca,'XScale','log')
set(gca,'YScale','log')
hold on; set(gcf,'Color','white'); box on

plot(x_vector,J,'Linewidth',1.5)

% Set the x-axis quantity
xlabel('[H_2SO_4] (cm^{-3})')
%xlabel('[NH_3] (cm^{-3})')
%xlabel('{\itT} (K)')
ylabel('{\itJ} (cm^{-3} s^{-1})')

