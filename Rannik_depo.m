function [dep_veloc] = Rannik_depo(T, dp, dens, U)
%	Semi-empirical deposition model presented by Rannik et al., JGR, 2003
%	* Note:  valid in the diameter range from around 10 to 500 nm 
%   * Note2: see the beginning of the code to see the values of fit
%            parameters
%
% INPUT:
% - T, temperature (in K) 
% - dens, particle density (kg/m3) 
% - dp, particle aerodynamic diameter (m) 
% - mean wind velocity (m/s) 
%
% OUTPUT: 
% dep_veloc - deposition velocity (m/s)
% 
% REFERENCE: 
% Rannik et al.: 
% "Interpretation of aerosol particle fluxes over a pine forest: Dry
% deposition and random errors", 
% J. Geophys. Res., 108, NO. D17, 4544, doi:10.1029/2003JD003542, 2003

d_min= 10e-9; d_max= 500e-9; % in nm

% if the diameter is <d_min (nm) , use the value @ d_min: 
dp= max(d_min, dp); 
% if the diameter is >d_max (nm), use the value @ d_max: 
dp= min(d_max, dp); 
       
% key fit parameters, based on the Hyytiälä measurements (see fig. 7 in the paper):
dexponent = 2.5;		% (from Rannik et al., fig. 7; TRFE < 50%)
collection = 2.9;		% collection efficiency factor (from fig. 7; TRFE < 50%)
dmin = 146d-9;			% diameter of minimum in collection efficiency (from fig. 7; TRFE < 50%)
%	dmin = 51.d-9
%	dexponent = 3.2

% other input parameters to the model
wind_speed = U;		    % median wind speed [m/s]
fric_veloc = 0.4;		% median friction velocity at Hyytiälä [m/s]
wind_canopy = 1.5;		% wind speed at canopy height [m/s]
gamma = 1.9;			% coefficient of exponential decrease of wind speed inside canopy

% some constants: 
K= 1.38D-23;
g= 9.81;
cv_to_cd = 1./3.; % the ratio of viscous to total drag (chosen according to Slinn, 1982) 
A2= 4.8*1d-4; 
radius= dp/2.0; % particle radius 

%	 ambient air properties
air_mean_path = (6.73d-8*T*(1.+110.4/T))/(296.*1.373);	 % mean free path
air_viscosity = (1.832d-5*406.4*T^1.5)/(5093*(T+110.4)); % viscosity
air_density = 1.2929*273.15/(T);		% [kg/m^3]     

%	Cunningham slip correction factor, eq. 6 in the paper 
corr_coeff = 1. + air_mean_path/(2.*radius)*(2.514 + 0.8*exp(-0.55*(2.*radius/air_mean_path)));
%	 particle settling velocity [m/s], eq. 5 in the paper
settling_veloc = (corr_coeff*dens*dp*dp*g)/(18.*air_viscosity);
%	Stokes number (not used below) 
ST= settling_veloc*fric_veloc/(g*A2);
%	 aerodynamic resistance, eq. 4 in the paper  [s/m]  
r_a = (wind_speed - wind_canopy)/(fric_veloc^2.0);

%	 collection efficiency of particles due to Brownian diffusion
diff_part = K*T*corr_coeff/(6.*3.41*air_viscosity*radius);	% eq. 12 in the paper
schmidt = air_viscosity/(air_density*diff_part); % Schmidt number, eq. 11
brownian = 1./3.*schmidt^(-2./3.);	% eq. 10							

% Calculate the prefactor in eq. 13 as a function of dmin 
coeff_dmin = 1. + air_mean_path/dmin*(2.514 + 0.8*exp(-0.55*dmin/air_mean_path)); % subtitute dmin into eq. 6
diff_dmin = K*T*coeff_dmin/(3.*3.41*air_viscosity*dmin); % substitute dmin into eq. 12 
schmidt_dmin = air_viscosity/(air_density*diff_dmin); % substitute dmin into eq. 11
deriv_coeff = -2.514*air_mean_path/(dmin^2) - (air_mean_path/dmin + 0.55)*0.8/dmin*exp(-0.55*dmin/air_mean_path); % deriv. of eq. 6
deriv_schmidt = schmidt_dmin*(1./dmin - 1./coeff_dmin*deriv_coeff); % eq. 15 
prefactor = cv_to_cd*(2./3.)*(schmidt_dmin^(-5./3.))*deriv_schmidt/(dexponent*dmin^(dexponent-1.)); % eq. 14
%	Finally, calculate empirical collection efficiency of particles, eq. 13 
empirical = 4.*prefactor*(radius^dexponent);
	
%	total collection efficiency of particles according to Rannik et al:
eps = collection*(brownian + empirical);	

%	 canopy resistance [s/m], eq.7 in the paper
eps = sqrt(eps);
r_c= wind_canopy/(fric_veloc^2*eps)*((1. + eps*tanh(gamma*eps))/(eps + tanh(gamma*eps)));
%	 deposition velocity [m/s]
dep_veloc = settling_veloc + 1./(r_a + r_c);	

end

