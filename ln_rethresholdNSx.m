function NEV = ln_rethresholdNSx(varargin)
% ln_rethresholdNSx_NEW
% opens headerless (*.dat) 30kHz datafile, filters it and finds threshold crossings
addpath('~/code/DataPreprocess')
addpath(genpath('~/matlab_functions/NPMK'))

[fileName pathName] = getFile;
kilorez = load([pathName,'rez.mat']);
NEV = openNEV();

%% Setting stuff up
newSpikes.TimeStamps = 0;
newSpikes.Electrode = 0;
newSpikes.Unit = 0;
newSpikes.Waveform = zeros(48,1);
firstTimeFlag = 1;
wavelengthLength = length(NEV.Data.Spikes.Waveform(:,1));

availableChannels = 32;

%% Highpass filter continuous data
hp_traces = batch_filter_trace([pathName,fileName],kilorez.rez.ops,'high',3,300);
trace_length = size(hp_traces,2);

%% Finding the spikes
for channelIDX = 1:availableChannels
    fprintf('Now processing channel %d\n',channelIDX)
    my_threshold = -4*median((abs(hp_traces(channelIDX,:))./0.6745));
    thresholdTimestamps = find(diff(hp_traces(channelIDX,:) < my_threshold) == 1);

    % find spikes triggered by the same event and take the first
    ovrlp_events = find(diff(thresholdTimestamps) < wavelengthLength);
    thresholdTimestamps(ovrlp_events) = [];
    
    inds = find((trace_length - thresholdTimestamps) < wavelengthLength);
    thresholdTimestamps(inds) = [];
    % clip if spike started before the recording 
    inds = find((thresholdTimestamps - 10) < 1);
    thresholdTimestamps(inds) = [];
    
    % while length(hp_traces(channelIDX, :)) - thresholdTimestamps(end) < 38
    %     thresholdTimestamps(end) = [];
    % end

    newSpikes.TimeStamps(end+1:end+length(thresholdTimestamps)) = thresholdTimestamps;
    newSpikes.Electrode(end+1:end+length(thresholdTimestamps)) = repmat(channelIDX, 1, length(thresholdTimestamps));
    newSpikes.Unit(end+1:end+length(thresholdTimestamps)) = zeros(1, length(thresholdTimestamps));
    if firstTimeFlag
        newSpikes.TimeStamps(1) = [];
        newSpikes.Electrode(1) = [];
        newSpikes.Unit(1) = [];
    end
    waveformIndices = bsxfun(@plus, thresholdTimestamps', -10:37)';
    waveformFlat = hp_traces(channelIDX, waveformIndices(:));
    newSpikes.Waveform(:,end+1:end+length(thresholdTimestamps)) = reshape(waveformFlat, 48, length(waveformFlat)/48);
    if firstTimeFlag
        newSpikes.Waveform(:,1) = [];
        firstTimeFlag = 0;
    end
    clear waveformFlat;
end

[~, timestampIndex] = sort(newSpikes.TimeStamps);

NEV.Data.Spikes.TimeStamp = newSpikes.TimeStamps(timestampIndex);
NEV.Data.Spikes.Electrode = newSpikes.Electrode(timestampIndex);
NEV.Data.Spikes.Unit = newSpikes.Unit(timestampIndex);
NEV.Data.Spikes.Waveform = newSpikes.Waveform(:,timestampIndex);

ln_saveNEV(NEV,[pathName,'thresholdedNEV.nev']);
