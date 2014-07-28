function [ out ] = StartConc( m,dp,d_limit,rho )
% StartConc Calculates concentration in the beginning of the simulation
% using the following parameters specifying the aerosol.
%
% m       = Total mass of particles released         (kg)
% dp      = Primary particle diameter                (m)
% d_limit = Size of the box of initial concentration (m)
% rho     = Particle density                         (kg/m^3)
%
% out     = Concentration                            (#/cm^3)
%
% Paxton Juuti & Joni Kalliokoski
% TTY 24.07.2014 
% 

m_p = rho*pi/6*dp^3;        % Single particle mass
N_p = m / m_p;              % Number of particles
out = N_p/(d_limit^3*1e6);  % Concentration in the initial dispersed volume

end
