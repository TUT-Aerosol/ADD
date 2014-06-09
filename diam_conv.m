function [d_aero] = diam_conv(T, d_ve, khi, rho)
% diam_conv(T, d_ve, khi, rho) - converts volume equivalent diameter to aerodynamic diameter
% INPUT: 
% T - temperature (K)
% d_ve - volume eqv. diameter (in nanometers, NOTE!) 
% khi - shape factor 
% rho - particle density (g/cm3)
% 
% OUTPUT: 
% d_aero (in nanometers)

options= optimset('Display','off'); % supress output from fsolve
amfp = (6.73d-8*T*(1.+110.4/T))/(296.*1.373);	% mean free path of air
rho_e= rho/ 1.0;                                % effective density (rho is in g/cm3)
cc = 1. + amfp/(d_ve)*(2.514 + 0.8*exp(-0.55*(2.*d_ve/amfp))); % Cunningham slip correction factor
apu= sqrt((rho_e*cc)/khi)*d_ve; 

% relation between d_ve and d_aero (Hinds, 1999):  
dv_da = @( x, apu, amfp) x*(1. + amfp/(x)*(2.514 + 0.8*exp(-0.55*(2.*x/amfp)))) - apu; 
x= fsolve( @(x) dv_da(x, apu, amfp), d_ve, options ); d_aero= x; % solve the equation

% for checking purposes 
% (d_aero*(1. + amfp/(d_aero)*(2.514 + 0.8*exp(-0.55*(2.*d_aero/amfp)))) - apu)

end

