function [fileDurations , fileNames] = kilo_split_ln(directory)

w_dir = pwd;
cd(directory);
D     = dir('*.ns5');
fileNames = cell(size(D));

for i = 1:length(D)
    fileNames{i} = D(i).name;
end

%file durations...
fileDurations = NaN(numel(fileNames) , 1);
for i = 1 : numel(fileNames)
    
    tmpData = openNSx([pwd,'/',fileNames{i}],'noread');
    fileDurations(i) = tmpData.MetaTags.DataDurationSec * tmpData.MetaTags.SamplingFreq;
    
end

cd(w_dir);