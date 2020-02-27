Ntot = 1000;
dp = logspace(0,3,100)*10^-3; % µm
rho= 1000;
[N,M] = respiratoryDeposition(Ntot,dp,rho);

figure;
subplot(211)
semilogx(dp*10^3,N)
ylabel('N (#/s)')

subplot(212)
loglog(dp*10^3,M)
xlabel('Hiukkaskoko (nm)')
ylabel('M (kg/s)')
