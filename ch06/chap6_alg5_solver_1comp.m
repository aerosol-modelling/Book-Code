% Use a solver - accurate and fast
[T,C] = ode15s(@(t,c) dgde_1comp(t,c,K,E,S,Q),[0,t_final],c_init);

% See that the concentrations are not negative - this might happen for very stiff systems
if min(min(C)) < -1e-50
    error('Negative concentrations')
else
    C(C<0) = 0.0;
end