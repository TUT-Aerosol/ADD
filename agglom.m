function [out] = agglom(out, agglo, plume, alfa, k, dt)
% AGGLOM calculates agglomeration rate following Lehtinen et al., JCI, 1996 
% (Current version based on Euler's method) 
% 
% INPUT (at time t): 
% * N_p - average number of primary particles in an aggregate
% * v_a - average volume of an agglomerate (m3)
% * phi - total aerosol volume concentration (m3 / cm3 air)
% * alfa - see eq. 4 in Lehtinen et al, JCI, 1996
% * k - Boltzmann's constant (J / K)
% * T - temperature (K)
% * rho - particle density (g/cm3)
% * D_f - fractal dimension 
% * dt - time step (in s)
% 
% OUTPUT (at time t+dt):
% * N_tot - total agglomerate number concentration (/cm3) 
% * N_p - average number of primary particles in an aggregate
% * v_a - average volume of an agglomerate (m3)

% "kinetic prefactor" in eq. 3 in Lehtinen et al.: 
c= 0.5*alfa*sqrt((6.*k*plume.T)/(1e3*agglo.rho))*((3./(4.*pi))^(1./6.));
% Right-hand side of eq. 3 in Lehtinen et al.: 
% (note conversion /cm3 -> /m3 for phi) 
delta_va= c* (out.v_a^(1./6.))*(1e6*out.phi)*(out.N_p^(2./agglo.D_f - 2./3.));

% update variables: 
out.v_a= out.v_a + delta_va*dt;
out.N_tot= out.phi / out.v_a; % coagulation conserves volume 
out.N_p=  out.v_a / out.v_p;
out.d_ve= ((6.*out.v_a)/pi)^(1./3.); % vol. equivalent diameter  

end

