function make_binary_ln(binaryFileName)

w_dir = pwd;
addpath('/home/lauri/matlab_functions/andrew_kilosort_concat/');
addpath(genpath('/home/lauri/matlab_functions/NPMK/NSx Utilities/'));
%merge all .ns5 files in a directory into a single binary file to pass to
%kilosort

%.ns5 files are here...
directory = uigetdir('' , '.NS5 Data Directory');

%write binary file to same directory...
cd(directory)

%Grab the ns5 file names...
[fileDurations , fileNames] = kilo_split_ln(directory);

%concatenate everything into a single data file...
for j = 1 : numel(fileNames)
    
    if j == 1
        fid = fopen(binaryFileName , 'w');
    else
        fid = fopen(binaryFileName , 'a');
    end

    tmpDat = openNSx([directory,'/' fileNames{j}]);
    fwrite(fid , tmpDat.Data , 'int16');
    fclose(fid);
    
end

save('concatDataInfo.mat', 'fileDurations', 'fileNames');
cd(w_dir);