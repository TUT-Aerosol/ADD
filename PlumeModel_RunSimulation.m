function [] = PlumeModel_RunSimulation(Filename)
%------------------------------------------------------------------------------
% FUNCTION FOR PLUME MODEL EXECUTION AND DATA SAVING
%
% First one must specify the desired parameter values in files _Agglos.mat and _Plumes.mat
%
% Next the user initiates the simulation loop which then 
% calculates the data and saves it.
%
% The user must take into account the data files size and set the
% appropriate number of parts for the data file so that all available RAM is 
% not exhausted. 
%
% FILES NEEDED:
% Agglos.mat.mat; 
% Plumes.mat
% The Plume model; All other files required by the plume model
% 
% This file uses parallel function, which requires the user to activate 
% parallel computing if using older MATLAB versions:  
% "matlabpool local" or "matlabpool open" in command view does this.
% This speeds up the simulation. Simulation time is the initial time 
% divided by the CPU core count.
%-----------------------------------------------------------------------------

load([Filename '_Plumes.mat'])
load([Filename '_Agglos.mat'])

Variable_count=length(Agglos); %Total different initial conditions

% SetN=Variable_count;
% Parts=(numel(dp)*numel(m)*numel(d_limit));
% Results{SetN}=[]; 
Results{Variable_count}=[];

% Simulation and Data saving

Parts=1; %Number of parts the data is split into

SetN=Variable_count/Parts;
tic %Total time of the simulation is displayed at the end to estimate the time taken by larger data sets
for ipart=1:Parts,
    parfor i=(ipart-1)*SetN+1:ipart*SetN;
        try
            [out]= agglo_disp_driv(Agglos{i}, Plumes{i});
            Results{i}=out;
            disp(['Loop number ' num2str(i) ' done']) %Displays the completed loops number for progress monitoring
        catch
            disp(['Loop number ' num2str(i) ' ERROR!!!!!!!!!!!!!']) 
            % Certain cases shown to produce errors, 
            % this prevents the simulation from ending due to single error
        end
    end
    save([Filename '_' num2str(ipart) '_Results.mat'], 'Results','-v7.3') %Save results
    clear Results; %Clear the results to free RAM
    Results{SetN*Parts}=[];
end
toc
disp(['*** Run for file: ' Filename ' DONE ***'])
end
