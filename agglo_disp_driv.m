function [out]= agglo_disp_driv(agglo, plume)
% agglo_disp_driv - Initializes and runs code simulating aerosol population
%                   undergoing agglomeration, dilution and deposition due to
%                   dispersion. Source is assumed to be at (0,0,0), and air 
%                   parcel is assumed to follow linear trajectory between 
%                   initial abd final coordinates. cartesian coordinates are
%                   aligned so that the x coordinate corresponds to wind
%                   direction.
%
% INPUT (contained in 'agglo' structure): 
% * N_0 - initial number concentration of primary particles(/cm3)
% * d_0 - diameter of the primary particles composing the agglomerates (m)
% * D_f - fractal dimension 
% * rho - particle density (g/cm3)
% * khi - particle shape factor
%
% INPUT VARIABLES (contained in 'plume' structure)
% * disp_scheme - dispersion scheme applied (choices: 'klug' and 'davidson')
% * depo_scheme - deposition scheme applied (choices: 'rannik')
% * U - wind velocity (m/s)
% * stab_class - stability class
% * N_reso - determines the time step (see below)
% * x_0, y_0, y_0 - coordinates for initial location
% * x_1, y_1, z_1 - coordinates for final location (trajectory is assumed
%   to follow linear path between these and initial coordinates)
% * T - temperature (K)
% * d_limit - neglect dilution, deposition and agglomeration at smaller 
%             distances than d_limit from the source (in m)
% 
% MAIN OUTPUT VARIABLES (contained in 'out' structure):
% * Ntot_ts - time series of number concentration (/cm3) 
% * va_ts - time series of average single agglomerate volume (m3) 
% * phi_ts - time series of total aerosol volume  (m3 / cm3 air)
% * Np_ts - time series of number of primary particles / agglomerate
% * x_ts, y_ts, z_ts, tc_ts - time series of spatial coordinates and time
%   variable (in m and s, respectively)
% * dist_ts - time series of distance from the source
% * sigx_ts, sigy_ts, sigz_ts - time series of dispersion parameters
% * Kl_ts - time series of 1st order loss parameter due to deposition and
%           dilution (in /s)
% * depvel_ts - time series of deposition velocity (in m/s) 
%
% USES THE FOLLOWING ROOUTINES: 
% * plume_conc.m - calculates the loss rate due to dilution and deposition
% * agglom.m - calculates the agglomeration rate and associated changes in
%   the aerosol properties
% * dp_klug.m / dp_davidson.m  - calculates dispersion parameters
% * Rannik_dep.m - calculates deposition rate (only choice at the moment,
%   4.4.2014)
% * diam_conv.m - diameter conversion, needed for deposition code  

% intervals of x, y and z coordinates corresponding to a single time step: 
dx= (plume.x_1 - plume.x_0); dy= (plume.y_1 - plume.y_0); dz= (plume.z_1 - plume.z_0); 
% total distance traversed
dist_tot= sqrt(dx^2 + dy^2 + dz^2); 
% time step (note: time step should be less than a second at the very least) 
out.dt= dist_tot/(plume.U*plume.N_reso); 

% Boltzmann's constant
k= 1.3806488e-23;

% time series of output variables:
out.x_ts= []; out.y_ts= []; out.z_ts= []; out.tc_ts= []; out.dist_ts= []; 
out.Ntot_ts= []; out.phi_ts=[]; out.va_ts= []; out.Np_ts= []; 
out.sigx_ts= []; out.sigy_ts= []; out.sigz_ts= []; out.Kl_ts=[]; out.depvel_ts= []; 
out.cx_ts= []; out.cy_ts= [];out.cz_ts= []; out.dcx_ts= []; out.dcy_ts= [];out.dcz_ts= [];
out.Ctest_ts= []; % for testing purposes

% calculate the alfa constant (eq. 4 in Lehtinen et al., JCI, 1996): 
alfa = 6.548 + 112.1*(agglo.D_f^(-7.883));

% set initial values

% initial total volume (phi_0) and volume of a primary particle (v_p)
out.phi0= agglo.N_0 * ((1./6.)*pi)*(agglo.d_0^3.0); % in m3 / cm3 air
out.v_p= out.phi0 / agglo.N_0; out.phi= out.phi0; % in m3
out.N_tot0 = agglo.N_0; out.N_tot= out.N_tot0; % initial particle concentration
out.v_a0 = out.v_p; out.v_a= out.v_a0; % initial average agglomerate volume
out.d_ve= ((6.*out.v_a)/pi)^(1./3.);
out.N_p = 1.0; % initially, only one primary particle / agglomerate
out.x= plume.x_0; out.y= plume.y_0; out.z= plume.z_0;  % initial coordinates
out.dist= sqrt(out.x^2 + out.y^2 + out.z^2); % distance from the source
out.U= plume.U; % wind velocity
out.tc= 0.0; % set current time 
out.disp_scheme_flag= true;

% ** main loop starts here **  
for i= 0: plume.N_reso
    % agglomeration: 
    if out.dist >= plume.d_limit
        out = agglom(out, agglo, plume, alfa, k, out.dt); 
    end
    
    % deposition:
    depo_scheme_flag= true;
    switch lower(plume.depo_scheme)
        case('rannik')
            % convert vol. eqv. to aerodynamic diameter: 
            out.d_a= 1e-9*diam_conv(plume.T, 1e9*out.d_ve, agglo.khi, agglo.rho); 
            % calculate deposition rate (in m/s): 
            out.dep_velo= Rannik_depo(plume.T, out.d_a, agglo.rho, plume.U);
        otherwise 
            depo_scheme_flag= false;  % unknown deposition scheme
            out.dep_velo= 0.0;
    end % end switch

    % calculate dispersion parameters: 
    switch lower(plume.disp_scheme)
            case('klug')
                dp= dp_klug(plume.stab_class, out.x); 
            case('davidson')
                dp = dp_davidson(plume.stab_class, out.x); 
            otherwise                 
                out.disp_scheme_flag= false;  % unknown dispersion scheme
                dp.stab_flag = false;
                dp.sig_x= 0.0; dp.sig_y= 0.0; dp.sig_z= 0.0;
                dp.dsig_x= 0.0; dp.dsig_y= 0.0; dp.dsig_z= 0.0;
    end % end switch
    
    % calculate loss rate due to dispersion and deposition: 
    if(dp.stab_flag == true && out.disp_scheme_flag == true && out.dist >= plume.d_limit) 
        dpout = plume_conc(out, dp);
    else
        dpout.K_loss= 0.0; dpout.c_x= 0.0; dpout.c_y= 0.0; dpout.c_z= 0.0; 
        dpout.dc_x= 0.0; dpout.dc_y= 0.0; dpout.dc_z= 0.0; dpout.test_conc= 0.0;
    end % endif
       
    % update volume and number concentration (Euler method):
    % out.N_tot= out.N_tot + out.N_tot*dpout.K_loss*out.dt;
    % out.phi= out.phi + out.phi*dpout.K_loss*out.dt;
    
    % update volume and number concentration (implicit Euler method):
    out.N_tot= out.N_tot/(1.0 - dpout.K_loss*out.dt);
    out.phi= out.phi/(1.0 - dpout.K_loss*out.dt);
   
% update time series of output variables: 
    out.x_ts= [out.x_ts out.x]; out.y_ts= [out.y_ts out.y];  out.z_ts= [out.z_ts out.z];  out.tc_ts= [out.tc_ts out.tc]; out.dist_ts= [out.dist_ts out.dist];     
    out.sigx_ts= [out.sigx_ts dp.sig_x]; out.sigy_ts= [out.sigy_ts dp.sig_y];  out.sigz_ts= [out.sigz_ts dp.sig_z];
    out.Ntot_ts= [out.Ntot_ts out.N_tot]; out.va_ts= [out.va_ts out.v_a]; out.phi_ts= [out.phi_ts out.phi]; out.Np_ts= [out.Np_ts out.N_p]; 
    out.Kl_ts= [out.Kl_ts abs(dpout.K_loss)]; out.depvel_ts= [out.depvel_ts out.dep_velo]; 
    out.cx_ts= [out.cx_ts dpout.c_x]; out.cy_ts= [out.cy_ts dpout.c_y]; out.cz_ts= [out.cz_ts dpout.c_z]; 
    out.dcx_ts= [out.dcx_ts dpout.dc_x]; out.dcy_ts= [out.dcy_ts dpout.dc_y]; out.dcz_ts= [out.dcz_ts dpout.dc_z]; 
    out.Ctest_ts= [out.Ctest_ts dpout.test_conc]; % for testing
    
% update consumed time (tc), local coordinates and distance from the source (dist):
    out.tc= out.tc + out.dt;
    out.x= out.x + dx/plume.N_reso; out.y= out.y + dy/plume.N_reso; out.z= out.z + dz/plume.N_reso;   
    out.dist= sqrt((out.x - plume.x_0)^2 + (out.y - plume.y_0)^2 + (out.z - plume.z_0)^2);     
end 
% ** main loop ends here **

if (out.disp_scheme_flag == false || dp.stab_flag == false)
    disp('** Unknown scheme for dispersion parameters or unknown stability class - dispersion and deposition were neglected **'); 
end 
if (depo_scheme_flag== false) 
    disp('*** Unknown deposition scheme - deposition was neglected ***');  
end 
    
% display output:
ADD_output(out, dpout); 

end
    






