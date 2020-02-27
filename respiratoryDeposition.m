function[NrespDep,MrespDep,LDSA] = respiratoryDeposition(Ntot, dp, rho)
% This function calculates human respiratory deposition of aerosol nanoparticles. Different respiratory
% regions concerned are head airways (ha), tracheo-bronchial (tb) and alveolar (al).
% Unimodal size distribution with size dp is assumed. Hinds suggest to use aerodynamic diameter for
% particles >0.5µm and physical/volume based diameter for particles <0.5µm.
% For more details see (Hinds, 1999) Aerosol Technology. ch. 11 pp. 233-245.
%
% INPUT:
% Ntot       total particle number concentration (#/cm3)
% dp         particle diameter (µm)
% rho        particle density (kg/m3)
%
% OUTPUT:
% NrespDep   total rate of deposited number of particles (#/s)
% MrespDep   total rate of deposited particle mass (kg/s)
% LDSA       lung deposited surface area per inhaled volume (µm2/cm3)
%
% Mikko Poikkimäki 4.11.2016

PLOTTING = 0; % plot DF:s on=1 off=0

% Volume inhaled by the exposed person per time unit (cm^3/s) (Hinds table 11.3)
% Female and male average of sitting, light and heavy exercise. 
Vm = (0.39+1.25+2.7+0.54+1.5+3)/6*10^6/3600;

% Inhalable fraction of particles (eq. 11.2)
IFi = 1-0.5.*(1-1./(1+0.00076.*dp.^2.8));
% Total deposition fraction (eq. 11.5)
DF = IFi.*(0.0587+0.911./(1+exp(4.77+1.485.*log(dp)))+0.943./(1+exp(0.508-2.58.*log(dp))));
% Tracheo-bronchial deposition fraction (eq. 11.3)
DFtb = 0.00352./dp.*(exp(-0.234.*(log(dp)+3.4).^2)+63.9.*exp(-0.819.*(log(dp)-1.61).^2));
% Alveolar deposition fraction (eq. 11.4)
DFal = 0.0155./dp.*(exp(-0.416.*(log(dp)+2.84).^2)+19.11.*exp(-0.482.*(log(dp)-1.362).^2));
% Head airways deposition fraction (eq. 11.1)
DFha = IFi.*(1./(1+exp(6.84+1.183.*log(dp)))+1./(1+exp(0.924-1.885.*log(dp))));

% Total rate of deposited number of particles per time unit (s)
NrespDep =  Ntot.*Vm.*DF;
% Total rate of deposited particle mass per time unit (s)
MrespDep = pi()./6.*rho.*(dp.*10^-6).^3.*NrespDep;
% Total LDSA per inhaled volume (µm2/cm3)
LDSA = pi().*dp.^2.*Ntot.*DF;

%TEST and PLOT DF:s
if PLOTTING == 1;
    %DFtest = DFha+DFtb+DFal;
    %semilogx(dp*10^3,DF,'k',dp*10^3,DFtest,'r',dp*10^3,DFha,'g',dp*10^3,DFtb,'b',dp*10^3,DFal,'m')
    semilogx(dp*10^3,DF,'k',dp*10^3,DFha,'k--',dp*10^3,DFtb,'k-.',dp*10^3,DFal,'k.','LineWidth',1.5)
    legend('Kokonaisdepositio','Pään hengitystiet','Henkitorvi ja keuhkoputkisto','Keuhkorakkulat') % finnish
    %legend('Total','Head airways','Tracheobronchial','Alveoli') % english
    xlabel('Hiukkaskoko (nm)')
    ylabel('Keuhkodeposoituva osuus, DF')
    set(gca,'Fontsize',12)
% 
%     % Text boxes
%     TextBox1 = uicontrol('style','text');
%     set(TextBox1,'String','Kokonais')
%     set(TextBox1,'Units','characters')
%     TextBox1Position = get(TextBox1,'Position');
%     TextBox1Position = [50 25 20 2];
%     set(TextBox1,'Position',TextBox1Position)
%     s = TextBox1.FontSize;
%     TextBox1.FontSize = 12;

end

end