function MakeParamFiles(trNum, slNum, thickness)

direcName = string('/home/klaferri/Desktop/Research/ThermalModel/Trough'+string(trNum)+'/Slice'+string(slNum)+'/Lag_'+string(thickness)+'mm/');
mkdir(direcName)
copyfile('/home/klaferri/Desktop/Research/ThermalModel/scripts/baseParam.txt', direcName)
cd(direcName)
movefile 'baseParam.txt' 'paramFile_slope.txt'
addpath('/home/klaferri/Desktop/Research/ThermalModel/scripts/')
C = readStruct('paramFile_slope.txt')
C.initialLagThickness = thickness/1000;
save("paramFile_slope.mat", "C");

end
