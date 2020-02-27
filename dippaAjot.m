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

tic
for i = 1:length(Filename)
    file = Filename{i};
    PlumeModel_RunSimulation(file);
end

disp('FINISHED!')
toc
timeElapsed = toc;
save('timeFilename','timeElapsed','Filename','-v7.3')

