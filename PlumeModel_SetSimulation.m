%------------------------------------------------------------------------------
% SCRIPT FOR PLUME MODEL INITIALIZATION, EXECUTION AND DATA SAVING
%
% First one must specify the desired parameter values, either as single
% values or vector of values. Every single possible combination of these
% initial values is prepared to file structures for the plume model to use.
%
% Next the user initiates the simulation loop which then 
% calculates the data and saves it.
%
% The user must take into account the data files size and set the
% appropriate number of parts for the data file so that all available RAM is 
% not exhausted. 
%
% FILES NEEDED:
% ADD_input.mat;   For use as the base for initialization
% SimQuad.m;       Specifies the simulation end points on a circular
%                  perimiter.
% StartConc.m;     Calculates the starting concentration from given particle
%                  morphology data and mass.
% The Plume model; All other files required by the plume model
% 
% This file uses parallel function, which requires the user to activate 
% parallel computing if using older MATLAB versions:  
% "matlabpool local" or "matlabpool open" in command view does this.
% This speeds up the simulation. Simulation time is the initial time 
% divided by the CPU core count.
% 
% Paxton Juuti & Joni Kalliokoski
% TTY 24.07.2014
% 
%-----------------------------------------------------------------------------

clear all;
load ADD_input.mat %Loads parameter matrices, which get replaced by the user set values below.

SAVE=1; %Data saving to files, 1=On, 0=Off.
Filename='140724_IniParam'; %Filename for the initial parameters

%Parameters initialization values
m           = [1];                          %Initial material mass (kg)
dp          = [6]*1e-9;                     %Initial particle size (m)
d_limit     = [100];                        %Initial ploom size(m)
rho         = 1.6*1e3;                      %Material density(kg/m^3)
d_f         = 3;                            %Fractal dimensionality for agglomerate
khi         = 1;                            %Dynamic shape factor

% U and stab_class are paired with each other and should be have the same
% length n which produces n number of initial conditions, not n*n.
U           = [10];                         %Wind speed (m/s) 
stab_class  = ['b'];                        %Stability class;   A:very unstable,B:unstable,C:slightly unstable,
                                            %               D:neutral,E:slightly stable,F:stable
                                            
Nl          = 91;                           %Number of simulation lines; see SimQuad.m for additional info
Dl          = 20000;                        %Maximum simulation distance (m); see SimQuad.m for additional info
[X,Y]       = SimQuad(Dl,Nl);               %Endpoints for simulation lines
Sim_N       = 50;                           %Number of lines simulated away from the wind direction: 
                                            %   read how many of the N1- lines are actually simulated  

agglo.khi   = khi;                          %Particle shape factor
plume.N_reso= 1000;                         %Number of steps per simulation; determines the step size with Dl and U
plume.disp_scheme = 'klug';                 %Dispersion scheme
plume.Ntot_limit = 1;                       %If the simulated total concentration N_tot (1/cm3) is lower than this limit, 
                                            %       the simulation is stopped to save simulation time
steps = plume.N_reso + 1;                   %Length of result vector needed for plotting

% Initializing the input structs for all possible combinations of input
% parameters.

i_num=1; %Running count from 1 to number of different cases.
for idp=1:numel(dp),                                    %Loop for Dp
    for im=1:numel(m),                                  %Loop for m
        for idlim=1:numel(d_limit),                     %Loop for d_limit
            for iusc=1:numel(U),                        %Loop for U
                for idf=1:numel(d_f),                   %Loop for d_f
                    for irho=1:numel(rho),              %Loop for rho
                        Dir=1;  %Angle reset to 0 
                        for ip=i_num:(i_num+Sim_N-1),   %Loop for angles
                            Plumes{ip}            = plume;
                            Plumes{ip}.x_1        = X(Dir);
                            Plumes{ip}.y_1        = Y(Dir);
                            Plumes{ip}.d_limit    = d_limit(idlim);
                            Plumes{ip}.stab_class = stab_class(iusc);
                            Plumes{ip}.U          = U(iusc);
                            Dir = Dir+1;
                            Agglos{ip}            = agglo;
                            
                            %Initial concentration (#/cm^3)
                            Agglos{ip}.N_0        = StartConc(m(im),dp(idp),d_limit(idlim),rho(irho)); 
                            Agglos{ip}.d_0        = dp(idp);
                            Agglos{ip}.D_f        = d_f(idf);
                            Agglos{ip}.rho        = rho(irho);
                        end
                        i_num = i_num+Sim_N;
                    end
                end
            end
        end
    end
end

Variable_count=iusc*idf*irho*idp*im*idlim*Sim_N; %Total different initial conditions

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

save([Filename '_Agglos.mat'], 'Agglos','-v7.3')
save([Filename '_Plumes.mat'], 'Plumes','-v7.3')
toc

%% Distance to certain concentration
% For map overlay
% Remember to load the data from where to draw the pictures !!!!!!!!
file = [Filename '_' num2str(ipart) '_Results.mat'];
load(file)

Contour_count=Variable_count/Sim_N; %Number of plume contours
Lines=4; %Number of contourlines
Distances=logspace(6,3,Lines); %Concentration values for contourlines

Dist=linspace(0,Dl,steps);   %Radial distance vector 
clear N_dist;                %Clears N_dist if it is calculated with different values
N_dist_ini(Lines,Sim_N)=0;   %Cell initialization
N_dist{Contour_count}=[];    %Cell matrix initialization
for coni=1:Contour_count;    %Loop that goes through all plumes
    N_dist{coni}=N_dist_ini; 
    for in=1:Lines,          %Loop for all concentrations defined in Distances
        ind=1;
        for is=1+((coni-1)*Sim_N):coni*Sim_N, 
            %Loop tests every step along a concentration line until certain
            %value is reached and then it stores the radial distance value
            %and then exits the loop
            for id=1:steps,  
                if Results{is}.Ntot_ts(id)<=Distances(in),
                    N_dist{coni}(in,ind)=Dist(id);
                    ind=ind+1;
                    break
                end
            end
        end
    end
end

thstep=90/(Nl-1); %Angle stepsize
Angles=([0 linspace(1,size(N_dist{1},2)-1,size(N_dist{1},2)-1)]).*thstep; %Angles for different plumelines

%X and Y coordinates for all points in all plumelines
for i2=1:numel(N_dist),
    Yl{i2}=N_dist{i2}.*0;
    Xl{i2}=N_dist{i2}.*0;
    Zl{i2}=N_dist{i2}.*0;
    for i1=1:Lines,
        Yl{i2}(i1,:)=sind(Angles).*N_dist{i2}(i1,:);
        Xl{i2}(i1,:)=cosd(Angles).*N_dist{i2}(i1,:);
        Zl{i2}(i1,:)=Distances(i1)*ones(1,size(N_dist{i2}(i1,:),2));
    end
    Yl{i2}=[-Yl{i2}(:,end:-1:1) Yl{i2}];
    Xl{i2}=[Xl{i2}(:,end:-1:1) Xl{i2}];
    Zl{i2}=[Zl{i2} Zl{i2}];
end

figure
for ip=1:numel(Yl),
    plot3(Xl{ip}',(Yl{ip}')+((ip-1)*400),(Zl{ip}'))
    hold on
    %plot3(Xl{ip}',(-Yl{ip}')+((ip-1)*10000),Zl{ip}')
end
% set(gca,'Zscale','log')
axis tight
view(2)
legend(num2str(Distances'))
xlabel ('Plume Distance (m)')
ylabel ('Plume Width (m)')

%% Map overlay
% Map overlay requires a map picture file from which the user must define
% the pixel count per kilometer. Few examples from a pictures scalebar:
% "821 - 1032 = 1km = 211 pixels from kartta.jpg"
% "6 - 116 = 2km -> 1km = 55 pixels from KarttaTRE.jpg"
% Ground zero is the origin for the plume and should be set after first
% plotting of the new image file as the x and y coordinates can be checked
% with the data cursor tool. Lastly the image requires the number of
% distance rings to be plotted for estimating the distances from the image.
% Each ring is 2km wider than the last to achieve kilometer wide sections.

WindDir=30; %'degrees from zero angle (wind from east to west)
Plumen=1; %The plume number to be plotted

figure(3)
clf
%Read the desired map file to overlay the plume into

%-------------------------- Example map 1
rgb=imread('kartta.jpg'); % Map file
GroundZero=[920 290];     % Coordinate shift of plume origin
Metre=211/1000;           % Scaling of plume to the scale of the map
Rings=2;                  % Number of distance rings
%--------------------------

%-------------------------- Example map 2
% rgb=imread('KarttaTRE.jpg'); 
% GroundZero=[856 475];
% Metre=55/1000;
% Rings=7;
%--------------------------

image(rgb)
hold on

PlumeX1=-Metre.*Xl{Plumen}; %Scaling the plume to match the scale of the map
PlumeY1=Metre.*Yl{Plumen};

PlumeXR1=PlumeX1.*cosd(-WindDir)-PlumeY1.*sind(-WindDir); %Rotating the plume to match the direction of the wind
PlumeYR1=PlumeX1.*sind(-WindDir)+PlumeY1.*cosd(-WindDir);

plot(GroundZero(1)+PlumeXR1',PlumeYR1'+GroundZero(2),'LineWidth',2) %Shifting the plumes origin

plot(GroundZero(1),GroundZero(2),'or','MarkerFaceColor',[1 0 0]) %Red dot for ground zero

D=Metre*1000*20; % 20km mark to use as the end point for scale markings

plot(GroundZero(1)+[0 D],GroundZero(2)+[0 0],':') %East-West and South-North crosshair
plot(GroundZero(1)+[-D 0],GroundZero(2)+[0 0],':')
plot(GroundZero(1)+[0 0],GroundZero(2)+[0 D],':')
plot(GroundZero(1)+[0 0],GroundZero(2)+[-D 0],':')

for i=2:2:2*Rings, %Distance rings around the ground zero
    Circle(GroundZero(1),GroundZero(2),i*1000*Metre);
    text(GroundZero(1)-i*1000*Metre,GroundZero(2),[num2str(i/2) 'km'],'FontSize',15);
end
text(GroundZero(1),GroundZero(2),' \leftarrow Wind','FontSize',18); %Distance labels

set(gca,'xticklabel',[]) %Removing ticks as they represent just the number of pixels in the map file
set(gca,'xtick',[])
set(gca,'yticklabel',[])
set(gca,'ytick',[])
ratio=daspect; %Set the image aspect ratio to match the image files
daspect(ratio)

h = text(GroundZero(1),GroundZero(2),' \leftarrow Wind','FontSize',18); %Label for wind direction
set(h, 'rotation', WindDir) %Rotate the label to match wind direction

legend([num2str(Distances') [' #/cm3';' #/cm3';' #/cm3';' #/cm3']]) %Concentration line legend



















