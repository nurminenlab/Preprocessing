% converts probe data to a KS-2 probe file
addpath('C:\Users\lonurmin\Desktop\code\npy-matlab')

date = '2025-04-03';
gate = '_g0';
task = 'Xori';
animal = 'MM004-Wolfjaw';

config_path = ['F:\localDATA\Electrophysiology\',animal,'\',date,'\',task,gate,'\',task,gate,'_imec0\'];

[filenames, paths] = uigetfile('*.npy','Pick probe info', 'MultiSelect','on');

chanMap0ind = readNPY([paths,filenames{1}]);
chanMap = readNPY([paths,filenames{2}]);
connected = readNPY([paths,filenames{3}]);
xcoords = readNPY([paths,filenames{4}]);
ycoords = readNPY([paths,filenames{5}]);
shankInd = ones(size(chanMap));

name = ['NP_',num2str(length(xcoords)),'_channels'];
save([config_path,name,'.mat'],'chanMap','chanMap0ind','connected','xcoords','ycoords','name');

fprintf([num2str(length(shankInd)),' channels were recorded from\n']);