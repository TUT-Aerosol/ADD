function [out]= agglo_disp_driv(agglo, plume)
% agglo_disp_driv - Initializes and runs code simulating aerosol population
%                   undergoing agglomeration, dilution and deposition due to
%                   dispersion. Source is assumed to be at (0,0,0), and air 
%                   parcel is assumed to follow linear trajectory between 
%                   initial and final coordinates. Cartesian coordinates are
%                   aligned so that the x coordinate corresponds to wind
%                   direction.
%
% Coded by Tatu Anttila, 2014
%   -additions and corrections: Mikko Poikkimäki 6.11.2016
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
% * BLH - boundary layer height (in m, optional) 
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
% time step (note: time step should be less than a second at the very
% least). Now resolution calculated by the preset timestep
out.dt = plume.dt;
plume.N_reso = round(dist_tot/(plume.U*out.dt)); %out.dt= dist_tot/(plume.U*plume.N_reso); 
out.N_reso = plume.N_reso;

% Boltzmann's constant
k= 1.3806488e-23;

% time series of output variables:
out.x_ts= zeros(1,plume.N_reso+1); out.y_ts= zeros(1,plume.N_reso+1); out.z_ts= zeros(1,plume.N_reso+1);  
out.dist_ts= zeros(1,plume.N_reso+1); out.tc_ts= zeros(1,plume.N_reso+1);
out.Ntot_ts= zeros(1,plume.N_reso+1); out.phi_ts=zeros(1,plume.N_reso+1); out.va_ts= zeros(1,plume.N_reso+1); out.Np_ts= zeros(1,plume.N_reso+1); 
out.sigx_ts= zeros(1,plume.N_reso+1); out.sigy_ts= zeros(1,plume.N_reso+1); out.sigz_ts= zeros(1,plume.N_reso+1); 
out.Kl_ts=zeros(1,plume.N_reso+1); out.depvel_ts= zeros(1,plume.N_reso+1); 
out.cx_ts= zeros(1,plume.N_reso+1); out.cy_ts= zeros(1,plume.N_reso+1);out.cz_ts= zeros(1,plume.N_reso+1); 
out.dcx_ts= zeros(1,plume.N_reso+1); out.dcy_ts= zeros(1,plume.N_reso+1);out.dcz_ts= zeros(1,plume.N_reso+1);
out.Ctest_ts= zeros(1,plume.N_reso+1);  % for testing purposes
out.normalization_ts = zeros(1,plume.N_reso+1); out.N_depEnv_ts = zeros(1,plume.N_reso+1);
out.M_depEnv_ts = zeros(1,plume.N_reso+1); out.N_Env_ts = zeros(1,plume.N_reso+1); out.M_Env_ts = zeros(1,plume.N_reso+1);
out.N_respDep_ts = zeros(1,plume.N_reso+1); out.M_respDep_ts = zeros(1,plume.N_reso+1); out.d_ve_ts = zeros(1,plume.N_reso+1);
%out.N_corr_ts{plume.N_reso+1} = [];

% calculate the alfa constant (eq. 4 in Lehtinen et al., JCI, 1996): 
% TODO: MOVE THIS TO agglom.m
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
out.disp_scheme_flag= true; out.BLH_flag= true; 

out.x_0= 0;

% ** main loop starts here **  
for i= 1: plume.N_reso+1
    %disp(['Before ' num2str(out.N_tot)])
    % agglomeration: 
    if out.dist >= plume.d_limit
        out = agglom(out, agglo, plume, alfa, k, out.dt); 
    end
    %disp(['After ' num2str(out.N_tot)])
    
    % deposition:
    depo_scheme_flag= true;
    switch lower(plume.depo_scheme)
        case('rannik')
            % convert vol. eqv. to aerodynamic diameter: 
            out.d_a= 1e-9*diam_conv(plume.T, 1e9*out.d_ve, agglo.khi, agglo.rho); 
            % calculate deposition & settling velocities (in m/s): 
            [out.dep_velo, out.set_velo]= Rannik_depo(plume.T, out.d_a, agglo.rho*10e3, plume.U);
        otherwise 
            depo_scheme_flag= false;  % unknown deposition scheme
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
    
    % check if the boundary layer height is specified, 
    % if so, cap the value of sigma_z to prevent further dispersion along z axis 
    % when the plume has reached the inversion layer
    if (~isfield(plume, 'BLH'))
        out.BLH_flag= false; dp.BLH= 1e99;
    else
        dp.BLH= plume.BLH; dp.sig_z= min(dp.sig_z, plume.BLH);
        if (out.x_0 == 0 && dp.sig_z == plume.BLH)
            out.x_0= out.x;
        end 
    end   
    
    % calculate loss rate due to dispersion and deposition:
    if(dp.stab_flag == true && out.disp_scheme_flag == true && out.dist >= plume.d_limit) 
        dpout = plume_conc(out, dp, 1);
    else
        dpout.K_loss= 0.0; dpout.c_x= 0.0; dpout.c_y= 0.0; dpout.c_z= 0.0; 
        dpout.dc_x= 0.0; dpout.dc_y= 0.0; dpout.dc_z= 0.0; dpout.test_conc= 0.0;
        if(out.dist <= plume.d_limit) 
            dpout.test_conc= 1.0;
        end
    end % end
       
    % update volume and number concentration due to dispersion and deposition(Euler method):
    % out.N_tot= out.N_tot + out.N_tot*dpout.K_loss*out.dt;
    % out.phi= out.phi + out.phi*dpout.K_loss*out.dt;
    
    % update volume and number concentration due to dispersion and deposition(implicit Euler method):
    out.N_tot= out.N_tot/(1.0 - dpout.K_loss*out.dt);
    out.phi= out.phi/(1.0 - dpout.K_loss*out.dt);
    
    %%% Respiratory and environmental deposition:
    % choose in which points (x,y,z) on the trajectory environmental and respiratory deposition are calculated
    percent = logspace(-2,2,50); % percent of the total distance [0.1, 1, 5, 10, 15, 20, 25, 50, 75, 100]; 
    out.index = round(plume.N_reso.*percent./100)+1;
    % go through all the chosen points and see if one of them is reached
    depo_calc_flag = false;
    for j = 1:length(out.index)
        if (i == out.index(j))
            depo_calc_flag = true;
        end
    end
    % calculate deposited amount of material if at the chosen point
    out.longer = 500; % simulate more steps in time for these observation points % TODO: move this to depositedAmount.m- function
    if (depo_calc_flag == true)
        out = depositedAmount(out, dp, agglo, plume, i, dx, dy, dz);
        %disp(['*** Deposition calculated at distance ' num2str(out.dist) ' m. ***'])
    else
        out.normalization = 0.0; out.N_corr = []; out.N_depEnv = 0.0;
        out.M_depEnv = 0.0; out.N_Env = 0.0; out.M_Env = 0.0; out.N_respDep = 0.0;
        out.M_respDep = 0.0;
    end
    %%%
    
% update time series of output variables: 
    out.x_ts(i)= out.x; out.y_ts(i)= out.y;  out.z_ts(i)= out.z;  out.tc_ts(i)= out.tc; out.dist_ts(i)= out.dist;     
    out.sigx_ts(i)= dp.sig_x; out.sigy_ts(i)= dp.sig_y;  out.sigz_ts(i)= dp.sig_z;
    out.Ntot_ts(i)= out.N_tot; out.va_ts(i)= out.v_a; out.phi_ts(i)= out.phi; out.Np_ts(i)= out.N_p; 
    out.Kl_ts(i)= abs(dpout.K_loss); out.depvel_ts(i)= out.dep_velo; 
    out.cx_ts(i)= dpout.c_x; out.cy_ts(i)= dpout.c_y; out.cz_ts(i)= dpout.c_z; 
    out.dcx_ts(i)= dpout.dc_x; out.dcy_ts(i)= dpout.dc_y; out.dcz_ts(i)= dpout.dc_z; out.BLH = dp.BLH;
    out.Ctest_ts(i)= dpout.test_conc; % for testing
    
    out.normalization_ts(i) = out.normalization; out.N_depEnv_ts(i) = out.N_depEnv;
    out.M_depEnv_ts(i) = out.M_depEnv; out.N_Env_ts(i) = out.N_Env; out.M_Env_ts(i) = out.M_Env; out.N_respDep_ts(i) = out.N_respDep;
    out.M_respDep_ts(i) = out.M_respDep; out.d_ve_ts(i) = out.d_ve;
    %out.N_corr_ts{i} = out.N_corr;
    
% update consumed time (tc), local coordinates and distance from the source (dist):
    out.x= out.x + dx/plume.N_reso; out.y= out.y + dy/plume.N_reso; out.z= out.z + dz/plume.N_reso;   
    out.dist= sqrt((out.x - plume.x_0)^2 + (out.y - plume.y_0)^2 + (out.z - plume.z_0)^2);     
    out.tc= out.tc + out.dt;

% stop the simulation if number concentration is below limit value Ntot_limit (#/cm3) to save simulation time
% simulate still until 2km
    if (out.N_tot < plume.Ntot_limit && out.dist >= 2000)
        break
    end
end 
% ** main loop ends here **

if (out.disp_scheme_flag == false || dp.stab_flag == false)
    disp('*** Unknown scheme for dispersion parameters or unknown stability class - dispersion and deposition were neglected ***'); 
end

if (depo_scheme_flag== false)
    disp('*** Unknown deposition scheme - deposition was neglected ***');  
end 

if (out.BLH_flag== false) 
    disp('*** Boundary layer height neglected ***');  
end
% display output:
% ADD_output(out, dpout); 

end
    






