function [out] = plume_conc(in, din)
% plume_conc (in, din) - Gaussian "puff" model, including deposition. 
% Source is assumed to be at the ground level, and gravitational settling
% is neglected. 
% For details, see:  
% * Llewelyn , R.P.: An analytical model for the transport, dispersion and elimination of air
%   pollutants emitted from a point source, Atmos. Environ., 17 (1983), pp. 249–256.
% * Stockie, J.M.: The Mathematics of Atmospheric Dispersion Modeling, 
%   SIAM Review 53 (2011), No. 2, pp. 349–372
%
% INPUT (contained in the structure 'in'): 
% * x,y,z,t - spatial coordinates (in m) & time (in s) 
% * U - wind velocity (m/s) 
% * dep_velo - deposition velocity (m/s) 
% INPUT (contained in the structure 'din'): 
% * sig_i (i= x,y,z) - dispersion parameters (in m) 
% * dsig_i (i= x, y, z) - derivates of dispersion parameters 
%                         with respect to x (m/m)
%
% OUTPUT (contained in the stucture 'out'): 
% * C - concentration at (x,y,z,t), 
% * d_C - time derivative of C
% * c_i (i= x,y,z) - C= Q*(c_x*c_y*c_z) (in /m), see below 
% * dc_i (i= x,y,z) - time derivatives of c_i 
% * K_loss - first order loss coefficient due to dilution + deposition 
%            (in /s, see below)

% convert dsig_i to derivatives with respect to time (in m/m -> m/s):
dsig_x= din.dsig_x * in.U; dsig_y= din.dsig_y * in.U; dsig_z= din.dsig_z * in.U;

% Here 0.5*d(sigma_i^2)/dt = (dsigma_i/dt)*sigma_i= K_i by definition 
% (K_i is the eddy diffusion coefficient, units m2/s)  
K_x= dsig_x*din.sig_x; K_y= dsig_y*din.sig_y; K_z= dsig_z*din.sig_z;

% prevent further dilution along z-axis if the plume has reached boundary layer height 
if (din.sig_z> din.BLH && in.x > in.x_0)
   din.sig_z= din.BLH; % cap sigma_z value
   K_z= 0.5*in.U*(din.sig_z^2/(in.x - in.x_0)); % use this formulation when sigma is constant 
end

gamma= -(in.dep_velo - in.set_velo)/K_z; % [gamma] = 1/m

% deposition term (see below): 
dpt1= gamma*exp(gamma*(in.z + 0.5*gamma*din.sig_z^2));
dpt2= erfc((in.z + gamma*din.sig_z^2)/din.sig_z); 
out.dep_rate= dpt1*dpt2;

% time derivative of the deposition rate: 
d_dep_rate= dpt1*dsig_z*(dpt2*(gamma^2*din.sig_z*dsig_z) - (2./(sqrt(pi))*exp(-((in.z + gamma*din.sig_z^2)/din.sig_z)^2)*(-in.z/din.sig_z^2 + gamma)));

% calculate concentration variables:  
out.c_x= (1./(sqrt(2.0*pi)*din.sig_x))*exp(-((in.x - in.U*in.tc)^2)/(2.0*din.sig_x^2));
out.c_y= (1./(sqrt(2.0*pi)*din.sig_y))*exp(-(in.y^2)/(2.0*din.sig_y^2));
cz_apu1= 2.0*exp(-(in.z^2)/(2.0*din.sig_z^2)); 
out.c_z= (1./(sqrt(2.0*pi)*din.sig_z))*(cz_apu1 + out.dep_rate);

% analytical concentration:
out.test_conc= out.c_x*out.c_y*out.c_z; 

% calculate time derivatives of concentration variables:  
out.dc_x= -(out.c_x/din.sig_x)*(dsig_x + dsig_x*(((in.x - in.U*in.tc)^2)/(din.sig_x^2)) + ((in.x - in.U*in.tc)*in.U)/din.sig_x); 
out.dc_y= -(out.c_y/din.sig_y)*dsig_y*(1.0 + in.y^2/(din.sig_y^2));
out.dc_z= (dsig_z/din.sig_z)*(-out.c_z - (in.z^2/(sqrt(2.0*pi)*din.sig_z^3))*cz_apu1) + d_dep_rate/(sqrt(2.0*pi)*din.sig_z); 

% first order loss coefficient (/s): 
out.K_loss= (out.dc_x/out.c_x + out.dc_y/out.c_y + out.dc_z/out.c_z);

end

