%% DIPPA plottailua, Mikko Poikkimäki
%% DONE General
load 07112016sensitivityOfBLH_1_Results.mat
load 07112016sensitivityOfBLH_Agglos.mat
load 07112016sensitivityOfBLH_Plumes.mat
i=1;
figure;
plot(Results{i}.tc_ts,Results{i}.dist_ts,'k.')
UU = Results{i}.dist_ts./Results{i}.tc_ts;
plot(UU)
Results{i}.U
Results{i}.index
Results{i}.N_tot0
Results{i}.dt

figure;
plot(Results{i}.normalization_ts)

figure;
i=4;
plot(Results{i}.dist_ts,'.')

%% DONE ADD basic
i = 1;
ADD_output(Results{i})

%% DONE Run time
% ei dippaan
figure;
for i = 1:length(Results)
    semilogy(i,Results{i}.dt,'o')
    hold on;
end
xlabel('Run index')
ylabel('dt')
%% DONE 1) Sensitivity for all parameters - LOAD files
Filename{1}  = '07112016sensitivityOfDlimit';
Filename{2}  = '07112016sensitivityOfBLH';
Filename{3}  = '07112016sensitivityOfM';
Filename{4}  = '07112016sensitivityOfDp';
Filename{5}  = '07112016sensitivityOfDf';
Filename{6}  = '07112016sensitivityOfRho';
Filename{7}  = '07112016sensitivityOfKhi';
Filename{8}  = '07112016sensitivityOfU';
Filename{9}  = '07112016sensitivityOfStabClass';
Filename{10} = '07112016sensitivityOfDispSchemeKlug';
Filename{11} = '07112016sensitivityOfDispSchemeDavidson';
Filename{12} = '07112016sensitivityOfT';

% Load files
R{length(Filename)-1} = [];
A{length(Filename)-1} = [];
P{length(Filename)-1} = [];

k=1; % from which file to end
for i = k:length(Filename)
    file = Filename{i};    
    if (i == 11)
        %Both DispSchemes to same cell
        r = load([file '_1_Results.mat']);
        a = load([file '_Agglos.mat']);
        p = load([file '_Plumes.mat']);
        for j = 1:length(r.Results)
            R{k-1}.Results{length(R{k-1}.Results)+1} =r.Results{j};
            A{k-1}.Agglos{length(A{k-1}.Agglos)+1} =a.Agglos{j};
            P{k-1}.Plumes{length(P{k-1}.Plumes)+1} =p.Plumes{j};
        end
        
    else
        R{k} = load([file '_1_Results.mat']) ;
        A{k} = load([file '_Agglos.mat']) ;
        P{k} = load([file '_Plumes.mat']) ;
        k = k+1;
    end
end

% k=8;
% for i = k:length(Filename)
%     disp(['k ' num2str(k)])
%     disp(['i ' num2str(i)])
%     k = k+1;
% end

%% DONE 2) Sensitivity for all parameters - PLOT results
% needs the cells R, A and P generated above, R = Outputs, A,P = Inputs 
clearvars -except R A P

% ii = 1;
% i=2;
% R{ii}.Results{i}.dt
% ADD_output(R{ii}.Results{i})

% Distance
observationPoint = 200; % 200, 1000, ...
observationLimit = 0.01;

% Plotting style
size = 10;
width = 1;
fontsize = 14;
markerEdgeColor =                   {'k',       'k',   'k',   'k',   'k',   'k',   'k',   'k', 'k',          'k',           'k'};
marker  =                           {'ks',      'kh',  'k^',  'k.',  'ko',  'k<',  'kd',  'k>','k*',         'kx',          'kv'};
varName = matlab.lang.makeValidName({'d_limit', 'BLH', 'N_0', 'd_0', 'D_f', 'rho', 'khi', 'U', 'stab_class', 'disp_scheme', 'T'});
Color = [0.7,0.7,0.7];
markerFaceColor = {[0.9,0.9,0.9], [0.1,0.1,0.1], [0.8,0.8,0.8], [0.7,0.7,0.7], [0.6,0.6,0.6], [0.5,0.5,0.5], [0.4,0.4,0.4],...
                    [0.3,0.3,0.3], Color, Color, [0.2,0.2,0.2]};
legendInfo = {'d_{limit} (m)', 'BLH (m)', 'N_0 (#/cm^3)', 'd_0 (m)', 'D_f', '\rho (g/cm^3)', '\chi', 'U (m/s)',...
              'Stabiilisuusluokka','Dispersioparametrisaatio', 'T (K)'};
legendInfo2 = {'d_{limit} (%)', 'BLH (%)', 'N_0 (%)', 'd_0 (%)', 'D_f (%)', '\rho (%)', '\chi (%)', 'U (%)',...
              'Stabiilisuusluokka','Dispersioparametrisaatio', 'T (%)'};

% Define reference point to calculate relative difference
ref_flag = false; % if true, calculate and plot relative values of y compared to default point
refx_flag = false; % if true, calculate and plot relative values also of x compared to default point
iRef = [9, 1, 5, 7, 5, 3, 1, 9, 7, 1, 7]; % reference point index for different parameters
iRef45 = iRef + 1; % reference point index for different parameters for 45 degrees from wind direction

individ_flag = false; % if true, plot all cases on individual figures

figure;
range = [1, 3:8, 11:length(R)]; % ALL
%range = 3:7; % AGGLOS
%range = [1,2,8,11]; % PLUMES
for ii = range % this loop goes through all different one-at-the-time cases , 
                                % skip 9 and 10 because they have string
                                % input value  and skip 2 because too high values
    if (individ_flag == true)
        figure;
    end
    for i = 1:2:length(R{ii}.Results) % in this loop other parameters remain unchanged but the one that was changed, choose values to wind dir.
        if (~isempty(R{ii}.Results{i})) % if there is data
            
            % Right x value (model input) to be plotted
            if (isfield(A{ii}.Agglos{i},varName{ii}))
                x = A{ii}.Agglos{i}.(varName{ii});
                xRef = A{ii}.Agglos{iRef(ii)}.(varName{ii});
            elseif (isfield(P{ii}.Plumes{i},varName{ii}))
                x = P{ii}.Plumes{i}.(varName{ii});
                xRef = P{ii}.Plumes{iRef(ii)}.(varName{ii});
            else
                disp('Right field in initial value cells (Agglos or Plumes) not found.')
            end
            
            iDist = find(and(R{ii}.Results{i}.dist_ts >=(observationPoint-observationLimit),R{ii}.Results{i}.dist_ts <=(observationPoint+observationLimit))); % index of the vector in distance 1000 m
            iDistRef = find(and(R{ii}.Results{iRef(ii)}.dist_ts >=(observationPoint-observationLimit),R{ii}.Results{iRef(ii)}.dist_ts <=(observationPoint+observationLimit))); % index of the reference vector in distance 1000 m

            subplot(221)
            y = R{ii}.Results{i}.N_Env_ts;
            yRef = R{ii}.Results{iRef(ii)}.N_Env_ts;
            y(y==0) = nan;
            if (ref_flag == true)
                if (refx_flag == true)
                    semilogx(abs((x-xRef)/xRef*100),(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                else
                    semilogx(x,(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                end
                ylabel('\DeltaN_{maaperä} (%)')
            else
                semilogx(x,y(iDist)*400,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                ylabel('Maaperän lukumääräpit. N (#/m^3)')
            end
            set(gca,'FontSize', fontsize)
            hold on;

            subplot(222)
            y = R{ii}.Results{i}.M_Env_ts;
            yRef = R{ii}.Results{iRef(ii)}.M_Env_ts;
            y(y==0) = nan;
            if (ref_flag == true)
                if (refx_flag == true)
                    loglog(abs((x-xRef)/xRef*100),abs((y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100),marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                else
                    loglog(x,abs((y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100),marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                end
                ylabel('|\DeltaM_{maaperä}| (%)')
            else
                loglog(x,y(iDist)*1e9*400,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                ylabel('Maaperän massapit. M (µg/m^3)')
            end
            set(gca,'FontSize', fontsize)
            hold on;  

            subplot(223)
            y = R{ii}.Results{i}.N_respDep_ts;
            yRef = R{ii}.Results{iRef(ii)}.N_respDep_ts;
            y(y==0) = nan;
            if (ref_flag == true)
                if (refx_flag == true)
                    semilogx(abs((x-xRef)/xRef*100),(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                    if (individ_flag == true)
                        xlabel(['|\Delta' legendInfo2{ii} '|'])
                    else
                        xlabel('|\DeltaSyöttöarvo| (%)')
                    end
                else
                    semilogx(x,(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                    if (individ_flag == true)
                        xlabel(legendInfo{ii})
                    else
                        xlabel('Syöttöarvo')
                    end
                end
                ylabel('\DeltaN_{keuhkodepositio} (%)')
            else
                semilogx(x,y(iDist),marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                ylabel('Keuhkodepositio N (#)')
                if (individ_flag == true)
                    xlabel(legendInfo{ii})
                else
                    xlabel('Syöttöarvo')
                end
            end    
            set(gca,'FontSize', fontsize)
            hold on;

            subplot(224)
            y = R{ii}.Results{i}.M_respDep_ts;
            yRef = R{ii}.Results{iRef(ii)}.M_respDep_ts;
            y(y==0) = nan;
            if (ref_flag == true)
                if (refx_flag == true)
                    loglog(abs((x-xRef)/xRef*100),abs((y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100),marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                    if (individ_flag == true)
                        xlabel(['|\Delta' legendInfo2{ii} '|'])
                    else
                        xlabel('|\DeltaSyöttöarvo| (%)')
                    end
                else
                    loglog(x,abs((y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100),marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                    if (individ_flag == true)
                        xlabel(legendInfo{ii})
                    else
                        xlabel('Syöttöarvo')
                    end
                end
                ylabel('|\DeltaM_{keuhkodepositio}| (%)')
            else
                loglog(x,y(iDist)*1e9,marker{ii},'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{ii},'MarkerFaceColor',markerFaceColor{ii})
                ylabel('Keuhkodepositio M (µg)')                
                if (individ_flag == true)
                    xlabel(legendInfo{ii})
                else
                    xlabel('Syöttöarvo')
                end
            end    
            set(gca,'FontSize', fontsize)
            hold on;            
        end
    end
end

% Custom legend
if (individ_flag == false)
    h = zeros(length(range), 1);
    for ileg = 1:length(range)
        hold on;    
        h(ileg) = plot(0,0,marker{range(ileg)}, 'visible', 'off', 'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor{range(ileg)},'MarkerFaceColor',markerFaceColor{range(ileg)});
    end
    subplot(221)
    legend(h, legendInfo{range},'Location','northoutside','Orientation','Horizontal','Fontsize',fontsize+2);
end

%% DONE 3) Sensitivity for parameters in bar plots: BLH, stab_class, disp_scheme
% needs the cells R, A and P generated above, R = Outputs, A,P = Inputs 
% uses also other input from 2)

ref_flag = false; % relative values
faceColor = 'k';
norm_flag = false; % normalized axis
save_flag = false; % save figures

% Distance
observationPoint = 1000; % 200, 1000, ...
observationLimit = 0.01;

legendInfo = {'d_{limit} (m)', 'BLH (m)', 'm (kg)', 'd_0 (nm)', 'D_f', '\rho (g/cm^3)', '\chi', 'U (m/s)',...
              'Stabiilisuusluokka','Dispersioparametrisaatio', 'T (K)'};

%range = [2,9,10];
range = 1:11;
%range=2;
for ii = range        
    k=1;
    for i = 1:2:length(R{ii}.Results) % in this loop other parameters remain unchanged but the one that was changed, choose values to wind dir.
        if (~isempty(R{ii}.Results{i})) % if there is data            
                        
            iDist = find(and(R{ii}.Results{i}.dist_ts >=(observationPoint-observationLimit),R{ii}.Results{i}.dist_ts <=(observationPoint+observationLimit))); % index of the vector in distance 1000 m
            iDistRef = find(and(R{ii}.Results{iRef(ii)}.dist_ts >=(observationPoint-observationLimit),R{ii}.Results{iRef(ii)}.dist_ts <=(observationPoint+observationLimit))); % index of the reference vector in distance 1000 m

            y1(k) = R{ii}.Results{i}.N_Env_ts(iDist);
            yRef1(k) = R{ii}.Results{iRef(ii)}.N_Env_ts(iDistRef);

            y2(k) = R{ii}.Results{i}.M_Env_ts(iDist);
            yRef2(k) = R{ii}.Results{iRef(ii)}.M_Env_ts(iDistRef);

            y3(k) = R{ii}.Results{i}.N_respDep_ts(iDist);
            yRef3(k) = R{ii}.Results{iRef(ii)}.N_respDep_ts(iDistRef);

            y4(k) = R{ii}.Results{i}.M_respDep_ts(iDist);
            yRef4(k) = R{ii}.Results{iRef(ii)}.M_respDep_ts(iDistRef);
            
            k=k+1;     
        end
    end
    %%% Save the maximum change from default for each parameter
    maxValue1(ii) = max(abs((y1-yRef1)./yRef1*100));
    maxValue2(ii) = max(abs((y2-yRef2)./yRef2*100));
    maxValue3(ii) = max(abs((y3-yRef3)./yRef3*100));
    maxValue4(ii) = max(abs((y4-yRef4)./yRef4*100));
    %%%
    %%% Save the average change from default for each parameter
    avgValue1(ii) = mean(abs((y1-yRef1)./yRef1*100));
    avgValue2(ii) = mean(abs((y2-yRef2)./yRef2*100));
    avgValue3(ii) = mean(abs((y3-yRef3)./yRef3*100));
    avgValue4(ii) = mean(abs((y4-yRef4)./yRef4*100));
    %%%
    %%% Save the all values change from default for each parameter
    allValue1{ii} = abs((y1-yRef1)./yRef1*100);
    allValue2{ii} = abs((y2-yRef2)./yRef2*100);
    allValue3{ii} = abs((y3-yRef3)./yRef3*100);
    allValue4{ii} = abs((y4-yRef4)./yRef4*100);
    %%%
    %%% Save the all values absolute
    absValue1{ii} = y1*400;
    absValue2{ii} = y2*400;
    absValue3{ii} = y3;
    absValue4{ii} = y4;
    %%%
    
    h = figure;
    subplot(221) %%%%%%%%%     
    
    if (ref_flag == true)
        bar((y1-yRef1)./yRef1*100,'FaceColor',faceColor)
        ylabel('\DeltaN_{maaperä} (%)')
    else
        bar(y1,'FaceColor',faceColor)                
        ylabel('Maaperän lukumääräpit. N (#/m^3)')
    end    
    if (ii == 1)
        set(gca,'XTickLabel',{0.01, 0.1, 1, 10, 100, 1000})
    elseif (ii == 2)
        set(gca,'XTickLabel',{'10^{99}', '10^{4}', '10^{3}', '10^{2}', '10^{1}'})
    elseif (ii == 3)
        set(gca,'XTickLabel',{0.1, 1 10, 100, 1000})
    elseif (ii == 4)
        set(gca,'XTickLabel',{1,10,50,100,200,300,400,500})
    elseif (ii == 5) 
        set(gca,'XTickLabel',{1.5, 2, 2.5, 3, 3.5, 4}) % 1 deleted
    elseif (ii == 6 )
        set(gca,'XTickLabel',{0.5, 1.0, 1.5, 2, 5, 10, 20, 30})
    elseif (ii == 7)
        set(gca,'XTickLabel',{1, 1.1, 1.3, 1.5, 1.8, 2.0})
    elseif (ii == 8)
        set(gca,'XTickLabel',{0.5, 1, 2, 5, 10, 20, 40}) % 3, 30 deleted 
    elseif (ii == 9)
        set(gca,'XTickLabel',{'a', 'b', 'c', 'd', 'e', 'f'})
    elseif (ii == 10)
        set(gca,'XTickLabel',{'Klug', 'Davidson'})
    elseif (ii == 11)    
        set(gca,'XTickLabel',{223,248,273,298,323})
    end
    if (norm_flag == true)
        ylim([-100 50])
    end 

    subplot(222) %%%%%%%%%
    
    if (ref_flag == true)
        bar((y2-yRef2)./yRef2*100,'FaceColor',faceColor)
        ylabel('\DeltaM_{maaperä} (%)')
    else
        bar(y2*1e9,'FaceColor',faceColor)
        ylabel('Maaperän massapit. M (µg/m^3)')
    end
    if (ii == 1)
        set(gca,'XTickLabel',{0.01, 0.1, 1, 10, 100, 1000})
    elseif (ii == 2)
        set(gca,'XTickLabel',{'10^{99}', '10^{4}', '10^{3}', '10^{2}', '10^{1}'})
    elseif (ii == 3)
        set(gca,'XTickLabel',{0.1, 1 10, 100, 1000})
    elseif (ii == 4)
        set(gca,'XTickLabel',{1,10,50,100,200,300,400,500})
    elseif (ii == 5) 
        set(gca,'XTickLabel',{1.5, 2, 2.5, 3, 3.5, 4}) % 1 deleted
    elseif (ii == 6 )
        set(gca,'XTickLabel',{0.5, 1.0, 1.5, 2, 5, 10, 20, 30})
    elseif (ii == 7)
        set(gca,'XTickLabel',{1, 1.1, 1.3, 1.5, 1.8, 2.0})
    elseif (ii == 8)
        set(gca,'XTickLabel',{0.5, 1, 2, 5, 10, 20, 40}) % 3, 30 deleted 
    elseif (ii == 9)
        set(gca,'XTickLabel',{'a', 'b', 'c', 'd', 'e', 'f'})
    elseif (ii == 10)
        set(gca,'XTickLabel',{'Klug', 'Davidson'})
    elseif (ii == 11)    
        set(gca,'XTickLabel',{223,248,273,298,323})
    end    
    if (norm_flag == true)
        ylim([-40 40])
    end

    subplot(223) %%%%%%%%
    
    if (ref_flag == true) 
        bar((y3-yRef3)./yRef3*100,'FaceColor',faceColor)
        ylabel('\DeltaN_{keuhkodepositio} (%)')
    else
        bar(y3,'FaceColor',faceColor)
        ylabel('Keuhkodepositio N (#)')
    end        
    xlabel(legendInfo{ii})
    if (ii == 1)
        set(gca,'XTickLabel',{0.01, 0.1, 1, 10, 100, 1000})
    elseif (ii == 2)
        set(gca,'XTickLabel',{'10^{99}', '10^{4}', '10^{3}', '10^{2}', '10^{1}'})
    elseif (ii == 3)
        set(gca,'XTickLabel',{0.1, 1 10, 100, 1000})
    elseif (ii == 4)
        set(gca,'XTickLabel',{1,10,50,100,200,300,400,500})
    elseif (ii == 5) 
        set(gca,'XTickLabel',{1.5, 2, 2.5, 3, 3.5, 4}) % 1 deleted
    elseif (ii == 6 )
        set(gca,'XTickLabel',{0.5, 1.0, 1.5, 2, 5, 10, 20, 30})
    elseif (ii == 7)
        set(gca,'XTickLabel',{1, 1.1, 1.3, 1.5, 1.8, 2.0})
    elseif (ii == 8)
        set(gca,'XTickLabel',{0.5, 1, 2, 5, 10, 20, 40}) % 3, 30 deleted 
    elseif (ii == 9)
        set(gca,'XTickLabel',{'a', 'b', 'c', 'd', 'e', 'f'})
    elseif (ii == 10)
        set(gca,'XTickLabel',{'Klug', 'Davidson'})
    elseif (ii == 11)    
        set(gca,'XTickLabel',{223,248,273,298,323})
    end    
    if (norm_flag == true)
        ylim([-100 50])
    end

    subplot(224) %%%%%%%%%
    
    if (ref_flag == true)
        bar((y4-yRef4)./yRef4*100,'FaceColor',faceColor)                                    
        ylabel('\DeltaM_{keuhkodepositio} (%)')
    else
        bar(y4*1e9,'FaceColor',faceColor)
        ylabel('Keuhkodepositio M (µg)')             
    end        
    xlabel(legendInfo{ii})    
    if (ii == 1)
        set(gca,'XTickLabel',{0.01, 0.1, 1, 10, 100, 1000})
    elseif (ii == 2)
        set(gca,'XTickLabel',{'10^{99}', '10^{4}', '10^{3}', '10^{2}', '10^{1}'})
    elseif (ii == 3)
        set(gca,'XTickLabel',{0.1, 1 10, 100, 1000})
    elseif (ii == 4)
        set(gca,'XTickLabel',{1,10,50,100,200,300,400,500})
    elseif (ii == 5) 
        set(gca,'XTickLabel',{1.5, 2, 2.5, 3, 3.5, 4}) % 1 deleted
    elseif (ii == 6 )
        set(gca,'XTickLabel',{0.5, 1.0, 1.5, 2, 5, 10, 20, 30})
    elseif (ii == 7)
        set(gca,'XTickLabel',{1, 1.1, 1.3, 1.5, 1.8, 2.0})
    elseif (ii == 8)
        set(gca,'XTickLabel',{0.5, 1, 2, 5, 10, 20, 40}) % 3, 30 deleted 
    elseif (ii == 9)
        set(gca,'XTickLabel',{'a', 'b', 'c', 'd', 'e', 'f'})
    elseif (ii == 10)
        set(gca,'XTickLabel',{'Klug', 'Davidson'})
    elseif (ii == 11)    
        set(gca,'XTickLabel',{223,248,273,298,323})
    end      
    if (norm_flag == true)
        ylim([-40 40])
    end
    
    if (save_flag == true)
        fileN = ['sensitivityBarplot' varName{ii}];
        savefig(h,fileN)
        saveas(h,[fileN '.jpg']) 
        saveas(h,[fileN '.eps']) 
    end
    
    clear y1 y2 y3 y4 yRef1 yRef2 yRef3 yRef4
end

% plot maximum values
legendInfo3 = {'d_{lim}', 'B_{lh}', 'M', 'd_0', 'D_f', '\rho', '\chi', 'U', 'Sta','Di', 'T'};
h = figure;
subplot(221)
bar(maxValue1)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaN_{maaperä}| (%)')

subplot(222)
bar(maxValue2)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaM_{maaperä}| (%)')
set(gca,'YScale','log')

subplot(223)
bar(maxValue3)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaN_{keuhkodepositio}| (%)')

subplot(224)
bar(maxValue4)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaM_{keuhkodepositio}| (%)')
set(gca,'YScale','log')

if (save_flag == true)
fileN = 'sensitivityBarplotMaximumValues';
savefig(h,fileN)
saveas(h,[fileN '.jpg']) 
saveas(h,[fileN '.eps']) 
end

% plot average values
legendInfo3 = {'d_{lim}', 'B_{lh}', 'M', 'd_0', 'D_f', '\rho', '\chi', 'U', 'Sta','Di', 'T'};
h = figure;
subplot(221)
bar(avgValue1)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaN_{maaperä}| (%)')

subplot(222)
bar(avgValue2)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaM_{maaperä}| (%)')
set(gca,'YScale','log')

subplot(223)
bar(avgValue3)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaN_{keuhkodepositio}| (%)')

subplot(224)
bar(avgValue4)
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaM_{keuhkodepositio}| (%)')
set(gca,'YScale','log')

if (save_flag == true)
fileN = 'sensitivityBarplotAverageValues';
savefig(h,fileN)
saveas(h,[fileN '.jpg']) 
saveas(h,[fileN '.eps'])
end

whisker = 10;
% BARPLOT
legendInfo3 = {'d_{lim}', 'B_{lh}', 'M', 'd_0', 'D_f', '\rho', '\chi', 'U', 'Sta','Di', 'T'};
h = figure;
subplot(221)
C = [allValue1{range}];
grp = [zeros(1,length(allValue1{1})),ones(1,length(allValue1{2})),2*ones(1,length(allValue1{3})),3*ones(1,length(allValue1{4})),...
    4*ones(1,length(allValue1{5})),5*ones(1,length(allValue1{6})),6*ones(1,length(allValue1{7})),7*ones(1,length(allValue1{8})),...
    8*ones(1,length(allValue1{9})),9*ones(1,length(allValue1{10})),10*ones(1,length(allValue1{11}))];
boxplot(C,grp,'Whisker',whisker)
%gca.XAxis.TickLabelInterpreter = 'latex';
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaN_{maaperä}| (%)', 'Interpreter','tex')
ylim([0, 200])
set(gca,'Fontsize',13)

subplot(222)
boxplot([allValue2{range}],grp,'Whisker',whisker)
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaM_{maaperä}| (%)', 'Interpreter','tex')
set(gca,'YScale','log')
ylim([1, 10^10])
set(gca,'Fontsize',13)

subplot(223)
boxplot([allValue3{range}],grp,'Whisker',whisker)
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaN_{keuhkodepositio}| (%)', 'Interpreter','tex')
ylim([0, 200])
set(gca,'Fontsize',13)

subplot(224)
boxplot([allValue4{range}],grp,'Whisker',whisker)
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('|\DeltaM_{keuhkodepositio}| (%)', 'Interpreter','tex')
set(gca,'YScale','log')
ylim([10^-5, 10^10])
set(gca,'Fontsize',13)

if (save_flag == true)
fileN = 'sensitivityBarplotAllValues200';
savefig(h,fileN)
saveas(h,[fileN '.jpg']) 
saveas(h,[fileN '.eps']) 
end

% BARPLOT Absolute values
legendInfo3 = {'d_{lim}', 'B_{lh}', 'M', 'd_0', 'D_f', '\rho', '\chi', 'U', 'Sta','Di', 'T'};
h = figure;
subplot(221)
C = [absValue1{range}];
grp = [zeros(1,length(absValue1{1})),ones(1,length(absValue1{2})),2*ones(1,length(absValue1{3})),3*ones(1,length(absValue1{4})),...
    4*ones(1,length(absValue1{5})),5*ones(1,length(absValue1{6})),6*ones(1,length(absValue1{7})),7*ones(1,length(absValue1{8})),...
    8*ones(1,length(absValue1{9})),9*ones(1,length(absValue1{10})),10*ones(1,length(absValue1{11}))];
boxplot(C,grp,'Whisker',whisker)
%gca.XAxis.TickLabelInterpreter = 'latex';
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('Maaperän lukumääräpit. N (#/m^3)', 'Interpreter','tex')
%ylim([0, 3e5])
ylim([0, 1e4])
set(gca,'Fontsize',12)

subplot(222)
boxplot([absValue2{range}]*1e9,grp,'Whisker',whisker)
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('Maaperän massapit. M (µg/m^3)', 'Interpreter','tex')
set(gca,'YScale','log')
ylim([10^-5, 10^5])
set(gca,'Fontsize',12)

subplot(223)
boxplot([absValue3{range}],grp,'Whisker',whisker)
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('Keuhkodepositio N (#)', 'Interpreter','tex')
%ylim([0, 4*10^9])
ylim([0, 10^8])
set(gca,'Fontsize',12)

subplot(224)
boxplot([absValue4{range}]*1e9,grp,'Whisker',whisker)
set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel',legendInfo3)
ylabel('Keuhkodepositio M (µg)', 'Interpreter','tex')
set(gca,'YScale','log')
ylim([10^-5, 10^10])
set(gca,'Fontsize',12)

if (save_flag == true)
fileN = 'sensitivityBarplotAllAbsValues1000';
savefig(h,fileN)
saveas(h,[fileN '.jpg']) 
saveas(h,[fileN '.eps']) 
end

%% DONE Deposition for a reference distance for different time steps
load_flag = false;
if load_flag == true
    load 06112016sensitivityOfDt_1_Results.mat
    s=load('06112016sensitivityOfDt5s_1_Results.mat');
    Results{length(Results)+1} =s.Results{1};
end

% Distance
observationPoint = 1000; % 200, 1000, ...
observationLimit = 0.01;

fontsize = 12;
markerEdgeColor = 'k';
markerFaceColor = [0.5, 0.5, 0.5];
size = 10;
width = 1;
marker = 'ks';
marker2 = 'ko';
bar_flag = false; % plot barplot if true, normal plot if false
iRef = 4; % reference point index (dt=0.04)
ref_flag = true; % if true, calculate and plot relative values compared to default point

figure;
for i = 1:length(Results)
    
    iDist = find(and(Results{i}.dist_ts >=(observationPoint-observationLimit),Results{i}.dist_ts <=(observationPoint+observationLimit))); % index of the vector in distance 1000 m
    iDistRef = find(and(Results{iRef}.dist_ts >=(observationPoint-observationLimit),Results{iRef}.dist_ts <=(observationPoint+observationLimit))); % index of the reference vector in distance 1000 m
       
    subplot(221)
    y = Results{i}.N_Env_ts;
    yRef = Results{iRef}.N_Env_ts;
    y(y==0) = nan;
    if bar_flag == true
        bar(Results{i}.dt,y(iDist))
    elseif (ref_flag == true)
        semilogx(Results{i}.dt,(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)    
        ylabel('\DeltaN_{maaperä} (%)')
    else 
        semilogx(Results{i}.dt,y(iDist),marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)
        ylabel('Maaperän lukumääräpit. N (#/m^3)')
    end
    set(gca,'FontSize', fontsize)
    hold on;

    subplot(222)
    y = Results{i}.M_Env_ts;
    yRef = Results{iRef}.M_Env_ts;
    y(y==0) = nan;
    if bar_flag == true
        bar(Results{i}.dt,y(iDist)*1e9)
        ylabel('Maaperän massapit. M (µg/m^3)')
    elseif (ref_flag == true)
        semilogx(Results{i}.dt,(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)    
        ylabel('\DeltaM_{maaperä} (%)')
    else
        semilogx(Results{i}.dt,y(iDist)*1e9,marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)
        ylabel('Maaperän massapit. M (µg/m^3)')
    end
    set(gca,'FontSize', fontsize)
    hold on;  

%      [hAx,hLine1,hLine2] = plotyy(Results{i}.dt,y(iDist)*1e9,Results{i}.dt,(y(iDist)-yRef(iDist))./yRef(iDist),'semilogx','semilogx');
%         hLine1.Marker = 's';
%         hLine2.Marker = 'o';
% %         hLine1.Color = 'k';
% %         hLine2.Color = 'k';
%         hold on;
%     end
%     ylabel(hAx(1),'Maaperän massapit. M (µg/m^3)')
%     ylabel(hAx(2),'Ero oletusarvoon (%)')
%     hold on;  

    subplot(223)
    y = Results{i}.N_respDep_ts;
    yRef = Results{iRef}.N_respDep_ts;
    y(y==0) = nan;
    if bar_flag == true
        bar(Results{i}.dt,y(iDist))
    elseif (ref_flag == true)
        semilogx(Results{i}.dt,(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)    
        ylabel('\DeltaN_{keuhkodepositio} (%)')
    else
        semilogx(Results{i}.dt,y(iDist),marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)
        ylabel('Keuhkodepositio N (#)')
    end
    xlabel('Aika-askel (s)')
    set(gca,'FontSize', fontsize)
    hold on;
    
    subplot(224)
    y = Results{i}.M_respDep_ts;
    yRef = Results{iRef}.M_respDep_ts;
    y(y==0) = nan;
    if bar_flag == true
        bar(Results{i}.dt,y(iDist)*1e9)
    elseif (ref_flag == true)
        semilogx(Results{i}.dt,(y(iDist)-yRef(iDistRef))./yRef(iDistRef)*100,marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)    
        ylabel('\DeltaM_{keuhkodepositio} (%)')
    else
        semilogx(Results{i}.dt,y(iDist)*1e9,marker,'MarkerSize',size,'LineWidth',width,'MarkerEdgeColor',markerEdgeColor,'MarkerFaceColor',markerFaceColor)
        ylabel('Keuhkodepositio M (µg)')
    end    
    xlabel('Aika-askel (s)')
    set(gca,'FontSize', fontsize)
    hold on;
    %legendInfo{k} = ['Aika-askel ' num2str(Results{i}.dt) ' s']; k = k+1;
end

%legend('Ympäristö N (#/m^3)','Ympäristö M (kg/m^3)','Keuhkodepositio N (#)','Keuhkodepositio M (kg)')

%legend(legendInfo)

%% DONE Deposition for all distances
% also for all different input variations
load_flag = false;
if load_flag == true
    load 06112016sensitivityOfDt_1_Results.mat
    s=load('06112016sensitivityOfDt5s_1_Results.mat');
    Results{length(Results)+1} =s.Results{1};
end

load 09112016defaultValuesPlumeMapDt04sMoreDepoValues_1_Results.mat

size = 10;
width = 1;
marker = 'k.';
k=1;
fontsize = 12;

figure;
for i = 1%:length(Results)
    subplot(221)
    y = Results{i}.N_Env_ts;
    y(y==0) = nan;
    loglog(Results{i}.dist_ts,y,marker,'MarkerSize',size,'LineWidth',width)
    ylabel('Maaperän lukumääräpit. N (#/m^3)')
    set(gca,'FontSize', fontsize)
    hold on;

    subplot(222)
    y = Results{i}.M_Env_ts;
    y(y==0) = nan;
    loglog(Results{i}.dist_ts,y*1e9,marker,'MarkerSize',size,'LineWidth',width)
    ylabel('Maaperän massapit. M (µg/m^3)')
    set(gca,'FontSize', fontsize)
    hold on;

    subplot(223)
    y = Results{i}.N_respDep_ts;
    y(y==0) = nan;
    loglog(Results{i}.dist_ts,y,marker,'MarkerSize',size,'LineWidth',width)
    ylabel('Keuhkodeposoitio N (#)')
    xlabel('Etäisyys lähteestä (m)')
    set(gca,'FontSize', fontsize)
    hold on;
    
    subplot(224)
    y = Results{i}.M_respDep_ts;
    y(y==0) = nan;
    loglog(Results{i}.dist_ts,y*1e9,marker,'MarkerSize',size,'LineWidth',width)
    ylabel('Keuhkodeposoitio M (µg)')
    xlabel('Etäisyys lähteestä (m)')
    set(gca,'FontSize', fontsize)
    hold on;
    legendInfo{k} = ['Aika-askel ' num2str(Results{i}.dt) ' s']; k = k+1;
end

%legend('Ympäristö N (#/m^3)','Ympäristö M (kg/m^3)','Keuhkodepositio N (#)','Keuhkodepositio M (kg)')
legend(legendInfo)

%% DONE Plot moving plume
%load 07112016sensitivityOfDlimit_1_Results.mat
i=7;
clear legendInfo
%legendInfo{length(Results{i}.index)} = [];
%plotStyle = {'b-','b--','bx','g-','g--','gx','m-','m--','mx','c-','c--','cx','r-','r--','rx','k-','k--','kx'};
plotStyle = {'k-','k--','kx'};
size = 15;
width = 1.2;

k2=1;
k=1;
figure;
hold on;
for j = Results{i}.index
    if(~isempty(Results{i}.N_corr_ts{j}))
        semilogy(Results{i}.time_disp_ts,Results{i}.N_corr_ts{j},plotStyle{1},'MarkerSize',size,'LineWidth',width); k2=k2+1;
        semilogy(Results{i}.time_disp_ts,Results{i}.N_corr_ts{j}/Results{i}.normalization_ts(j),plotStyle{2},'MarkerSize',size,'LineWidth',width); k2=k2+1;
        semilogy(Results{i}.tc_ts(j),Results{i}.Ntot_ts(j),plotStyle{3},'MarkerSize',size,'LineWidth',width); k2=k2+1;
%         legendInfo{k} = ['C_{corr} (' num2str(Results{i}.dist_ts(j)) ' m)']; k = k+1;
%         legendInfo{k} = ['C_{gaussian}' '']; k = k+1;
%         legendInfo{k} = ['Referenssipiste C_{g,d}' '']; k = k+1;
    end
end

xlabel('Aika (s)')
ylabel('Lukumääräpitoisuus N (#/cm^3)')
axis([100 1200 1 10^6])
set(gca, 'YScale', 'log');
set(gca,'FontSize', 12)

%legend(legendInfo)
legend('C_{corr}', 'C_{gaussian}', 'Referenssipiste C_{g,d}')

% % Text boxes
TextBox1 = uicontrol('style','text');
set(TextBox1,'String','1 km')
set(TextBox1,'Units','characters')
%TextBox1Position = get(TextBox1,'Position');
%TextBox1Position = [15 20 15 1.7];
%set(TextBox1,'Position',TextBox1Position)
TextBox1.FontSize = 12;
TextBox1.BackgroundColor = 'w';

TextBox1 = uicontrol('style','text');
set(TextBox1,'String','2 km')
set(TextBox1,'Units','characters')
TextBox1.FontSize = 12;
TextBox1.BackgroundColor = 'w';

TextBox1 = uicontrol('style','text');
set(TextBox1,'String','3 km')
set(TextBox1,'Units','characters')
TextBox1.FontSize = 12;
TextBox1.BackgroundColor = 'w';

TextBox1 = uicontrol('style','text');
set(TextBox1,'String','4 km')
set(TextBox1,'Units','characters')
TextBox1.FontSize = 12;
TextBox1.BackgroundColor = 'w';

TextBox1 = uicontrol('style','text');
set(TextBox1,'String','5 km')
set(TextBox1,'Units','characters')
TextBox1.FontSize = 12;
TextBox1.BackgroundColor = 'w';

%semilogy(Results{i}.tc_ts(j),Results{i}.Ctest_ts(j)*Results{i}.N_tot0,'k*')

%% DONE Test analytic vs iterative conc
i=3;
figure;
loglog(Results{i}.dist_ts,Results{i}.Ctest_ts*Results{i}.N_tot0,'m.-')
%loglog(Results{i}.dist_ts,Results{i}.Ctest_ts*Results{i}.N_tot0/2.433e5*5.526e9,'k.-')
hold on;
semilogy(Results{i}.dist_ts,Results{i}.Ntot_ts,'k.-')
legend('Analytic with dispersion and deposition', 'Iterative with additional agglomeration') %'Analytic corrected',
xlabel('Etäisyys lähteestä (m)')
ylabel('Lukumääräpitoisuus (#/cm^3)')

%% DONE Test analytic vs iterative conc vs without agglomeration iterative
% DIPPAAN kuva!!!
i=1;
figure;

load 06112016sensitivityWithAgglomerationNreso100000dlimit001_1_Results.mat
loglog(Results{i}.dist_ts,Results{i}.Ntot_ts,'k--')
hold on;

load 06112016sensitivityWithAgglomerationNreso100000dlimit001WithSettlingVelocity_1_Results.mat
%loglog(Results{i}.dist_ts,Results{i}.Ntot_ts,'k-.')

loglog(Results{i}.dist_ts,Results{i}.Ctest_ts*Results{i}.N_tot0,'k.')

load 06112016sensitivityWithoutAgglomerationNreso100000dlimit001_1_Results.mat
semilogy(Results{i}.dist_ts,Results{i}.Ntot_ts,'k:')

loglog(Results{i}.dist_ts,Results{i}.Ctest_ts*Results{i}.N_tot0,'k')
%%%%loglog(Results{i}.dist_ts,Results{i}.Ctest_ts*Results{i}.N_tot0/2.433e5*5.526e9,'k.-')

for j = Results{i}.index
    if(~isempty(Results{i}.N_corr_ts{j}))
        semilogy(Results{i}.time_disp_ts(j),Results{i}.N_corr_ts{j}(j)/Results{i}.normalization_ts(j),'kx');
    end
end

% legend('Iterative with additional agglomeration','Iterative with additional agglomeration + settling velocity',...
%     'Analytic with dispersion and deposition + settling velocity',...
%     'Iterative with dispersion and deposition','Analytic with dispersion and deposition',...
%     'Analytic with only dispersion') %'Analytic corrected', ENGLISH

legend('1) Iteratiivinen:  agglomeraatiolla', '2) Analyyttinen: disp., dep. ja aset. nop.',...
    '3) Iteratiivinen:  disp. ja dep.','4) Analyyttinen: disp. ja dep.',...
    '5) Analyyttinen: vain disp.','Location','best') %'Analytic corrected', FINNISH

xlabel('Etäisyys lähteestä (m)')
ylabel('Lukumääräpitoisuus (#/cm^3)')
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca,'FontSize', 12)
