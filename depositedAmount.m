function [ out ] = depositedAmount( out, dp, agglo, plume, i, dx, dy, dz )
% This functions calculates deposited particle number and mass to human
% lung (respiratory) and to environment (soil). 
%
% Mikko Poikkimäki 5.11.2016

% number concentration N(x,y,z,t) without agglomeration or deposition (pure Gaussian puff model) for t= t_0 -> t_end
% REMARK: in here the time_disp runs through the whole simulation time
out.time_disp_ts = 0:out.dt:out.dt*(plume.N_reso+out.longer); % simulate time from the beginning of the dispersion event to the end for a single point (x,y,z)
N_disp_ts = zeros(1,plume.N_reso+1+out.longer);
N_disp = out.N_tot0;
for i2 = 1: plume.N_reso+1+out.longer
    %update time
    out.time_disp = out.time_disp_ts(i2);

    % calculate conc at (x,y,x,t) only dispersion
    % if distance is lower than d_limit then conc = initial conc until the
    % plume has blown away
    if (dp.stab_flag == true && out.disp_scheme_flag == true && out.dist <= plume.d_limit && out.time_disp <= plume.d_limit/plume.U)
        N_disp = out.N_tot0;
    elseif(dp.stab_flag == true && out.disp_scheme_flag == true) % && out.dist >= plume.d_limit) 
        dispOnly = plume_conc(out,dp,0);
        N_disp = out.N_tot0*dispOnly.test_conc; % analytical concentration
    else
        N_disp = 0; 
    end
        
    % update concentration
    N_disp_ts(i2) = N_disp; % concentration in the observation point x,y,z through whole simulation time 
    
end % end for

% Normalization coeff, difference between calculated number concentration with and without agglomeration and deposition
out.normalization = out.N_tot/N_disp_ts(i); % with depos. per without
out.N_corr = N_disp_ts*out.normalization; % corrected number concentration

% calculate deposited amount of particles per area (#/m2) on the ground
out.N_depEnv = trapz(out.time_disp_ts,out.N_corr*out.dep_velo); % deposited number per m2
out.M_depEnv = out.N_depEnv*out.v_a*agglo.rho*10e3; % deposited mass kg/m2

%assume that the mass/number deposited on the surface is mixed within depth (mixDepth) and distributed uniformly
out.mixDepth = 0.05; % m, value for urban and natural soil as used by Mueller and Nowack (2008) Env. Sci. Tech. 42,4447-4453 

% calculate environmental concentrations in soil
out.N_Env = out.N_depEnv/out.mixDepth; % number concentration #/m3 in soil
out.M_Env = out.M_depEnv/out.mixDepth; % mass concentration kg/m3 in soil

% calculate deposited amount of particles in human lung = respiratory deposition
[N_respDepRate,M_respDepRate] = respiratoryDeposition(out.N_corr,out.d_ve*10^6,agglo.rho*10e3); % Respiratory deposition rates #/s and kg/s
out.N_respDep = trapz(out.time_disp_ts,N_respDepRate); % deposited number to lungs
out.M_respDep = trapz(out.time_disp_ts,M_respDepRate); % deposited mass (kg) to lungs

end

%%%% TESTING

% % from now on out.x,y,z means distance between the plume and the observation point x,y,z
% % save the coordinates of the observation point
% x_disp = out.x; 
% y_disp = out.y;
% z_disp = out.z;

%     % calculate dispersion parameters: 
%     switch lower(plume.disp_scheme)
%         case('klug')
%             dp= dp_klug(plume.stab_class, out.x); 
%         case('davidson')
%             dp = dp_davidson(plume.stab_class, out.x); 
%         otherwise                 
%             out.disp_scheme_flag= false;  % unknown dispersion scheme
%             dp.stab_flag = false;
%             dp.sig_x= 0.0; dp.sig_y= 0.0; dp.sig_z= 0.0;
%             dp.dsig_x= 0.0; dp.dsig_y= 0.0; dp.dsig_z= 0.0;
%     end % end switch
%         
%     % check if the boundary layer height is specified, 
%     % if so, cap the value of sigma_z to prevent further dispersion along z axis 
%     % when the plume has reached the inversion layer
%     if (~isfield(plume, 'BLH'))
%         out.BLH_flag= false; dp.BLH= 1e99;
%     else
%         dp.BLH= plume.BLH; dp.sig_z= min(dp.sig_z, plume.BLH);
%         if (out.x_0 == 0 && dp.sig_z == plume.BLH)
%             out.x_0= out.x;
%         end 
%     end   

%     N_disp = N_disp/(1.0 - dispOnly.K_loss*out.dt);

%     % update the distance between the observation point and the plume
%     out.x = abs(out.x-dx/plume.N_reso);
%     out.y = abs(out.y-dy/plume.N_reso);
%     out.z = abs(out.z-dz/plume.N_reso);
    
% % Replace to original values (coordinates of this observation point)
% out.x = x_disp; 
% out.y = y_disp;
% out.z = z_disp;

