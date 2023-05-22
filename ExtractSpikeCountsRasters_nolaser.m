function ExtractSpikeCountsRasters_nolaser()
addpath(genpath('/home/lauri/matlab_functions/NPMK'));
addpath(genpath('/home/lauri/code/npy-matlab'));
addpath(genpath('/home/lauri/code/spikes'));
addpath('/home/lauri/matlab_functions/statisticals');

plot_all = false
%close all
errb = 1;
clc;

base_dir    = '/opt3/';
animal      = 'MK374/';
penetration = 'P2/';

% anal params
% in s
anal_offset   = 0.05;
anal_duration = 0.4;
base_duration = 0.4;
base_onset    = 0.1;
base_offset   = base_onset + anal_duration;
% in ms
raster_trail = 100;

% open NEV and get spikes
tmp_fld = [base_dir,animal,penetration];
sDir = uigetdir(tmp_fld);
sDir = [sDir, '/'];
NEV = openNEV([sDir,'thresholdedNEV.nev-rthr.nev'])
NEVfreq  = NEV.MetaTags.TimeRes;
spt = NEV.Data.Spikes.TimeStamp;
clu = NEV.Data.Spikes.Electrode;

stparf = dir([sDir,'stimParams*.mat']);
load([stparf(1).folder,'/',stparf(1).name])

% get pointer to dat-file
datf = dir(fullfile(sDir,'*-001.ns5.dat'));

% open continuous TTL channel record
a(1) = dir([sDir,'20*.ns4']);
DataStruct  = openNSxNew([a(1).folder,'/',a(1).name],'read','c:1');
Fs          = DataStruct.MetaTags.SamplingFreq;
% generate pulse sample times from continuous trigger channel
pulseTimes = PulseCounter(DataStruct.Data,20000);
pulseTimes = round(pulseTimes*(NEVfreq/Fs));
% because Maryam records one pulse at the beginning and one in the
% end of a block of stimuli
pulseTimes = pulseTimes(1:2:end) + (stimParams.off_time_ms/1000) * NEVfreq;

% get start time for each trial
stimTimes = NaN * ones(size(stimParams.stimOrder));
for j = 1:length(pulseTimes)
    for i = 1:size(stimTimes,1)
        stimTimes(i,j) = pulseTimes(j) + (i-1) * (((stimParams.off_time_ms + stimParams.on_time_ms)/1000) * NEVfreq);
    end
end

%% collect stimulus info
stimOnCount    = round((stimParams.on_time_ms/1000)*Fs);
stimOffCount   = round((stimParams.off_time_ms/1000)*Fs);
stimTotalCount = stimOnCount+stimOffCount;
stimMatrix     = stimParams.stimOrder;

nStims    = size(stimParams.stimOrder,1);
nTrials   = size(stimParams.stimOrder,2);

% trials and stim params to NEV time
stimTimesNEV  = round(stimTimes);
stimLengthNEV = round(stimOnCount*(NEVfreq/Fs));
stimOffNEV    = round(stimOffCount*(NEVfreq/Fs));

% Use time after last stimulus for blank response
blnkTimes = stimTimesNEV(end,:) + stimLengthNEV + 0.1*30e3;

N_units       = length(unique(clu));
units         = unique(clu);

raster_trail = 100;
raster_befre = 400;    

spikeCnts_NoL_tmp = NaN*ones(N_units, nStims, nTrials);
baseLine_tmp      = NaN*ones(N_units, nStims, nTrials);
spkraster_NoL_tmp = NaN*ones(N_units, nStims, nTrials, length(-raster_befre:1:stimParams.on_time_ms+raster_trail)-1);

for u = 1:N_units
    spikeTimeMat = double(spt(units(u) == clu));
    
    % collate spike-counts depending on weather the laser was on or off
    for f = 1:nStims
        idx = find(stimMatrix==f);        
        for g = 1:length(idx)

            % raster
            spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) - base_duration*30e3) & (spikeTimeMat<stimTimesNEV(idx(g)) + stimLengthNEV + base_duration*30e3));
            spk_tmp = (double(spikeTimeMat(spk_tmp) - stimTimesNEV(idx(g))))./30e3*1000;
            spkraster_NoL_tmp(u,f,g,:) = histcounts(spk_tmp,[-raster_befre:1:stimParams.on_time_ms+raster_trail]);
            % count
            spikeCnts_NoL_tmp(u,f,g) = sum(spikeTimeMat> (stimTimesNEV(idx(g)) + anal_offset*30e3) & spikeTimeMat<stimTimesNEV(idx(g)) + (anal_offset+anal_duration)*30e3);
        
            baseLine_tmp(u,f,g) = sum(spikeTimeMat > stimTimesNEV(idx(g))+stimLengthNEV+base_onset*30e3 & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV+base_offset*30e3);
        end
    end

    fprintf('Now fitting unit %d\n', u)
    if u == 1
        f = fit_sizetuning(mean(spikeCnts_NoL_tmp(u,:,:),3),unique(stimParams.VRF_Diam)');
        fit_params     = f.final_params;
        field_params   = f.field_params;
    else
        f = fit_sizetuning(mean(spikeCnts_NoL_tmp(u,:,:),3),unique(stimParams.VRF_Diam)');
        fit_params     = [fit_params; f.final_params];
        field_params   = [field_params; f.field_params];
    end

end

% for compatibility with downstream pytable
mn_spikeWidths = NaN*ones(size(units));
SNR       = NaN*ones(size(units));
waveForms = NaN*ones(size(units,1),10e3,82);
mxmp_cnt  = NaN*ones(size(units));
mn_spikeDepths = NaN*ones(size(units));

diams = unique(stimParams.VRF_Diam)';

% due to script history
spikeCnts_NoL = spikeCnts_NoL_tmp;
baseLine      = baseLine_tmp;
spkraster_NoL = spkraster_NoL_tmp;
contrast      = stimParams.Contrast*ones(size(spikeCnts_NoL_tmp,3),1);

writeNPY(spikeCnts_NoL, [sDir,'Kilosorted_spkCnts_NoL_400msMUA.npy']);
writeNPY(baseLine, [sDir,'Kilosorted_baseLine_400msMUA.npy']);
writeNPY(spkraster_NoL, [sDir,'Kilosorted_spkraster_NoL_400msMUA.npy']);
writeNPY(fit_params, [sDir,'Kilosorted_fitParams_NoL_400msMUA.npy']);
writeNPY(field_params, [sDir,'Kilosorted_fieldParams_NoL_400msMUA.npy']);
writeNPY(mn_spikeWidths, [sDir,'Kilosorted_spikeWidths_400msMUA.npy']);
% valid depth because batch_filter_trace.m reorganizes electrodes
writeNPY((1:size(spikeCnts_NoL,1))/10, [sDir,'Kilosorted_spikeDepths_400msMUA.npy']);
writeNPY(SNR, [sDir,'Kilosorted_SNR_400msMUA.npy']);
writeNPY(waveForms, [sDir,'Kilosorted_waveForms_400msMUA.npy']);
writeNPY(mxmp_cnt,[sDir,'Kilosorted_maxamp_contact_400msMUA.npy']);
writeNPY(diams, [sDir,'Kilosorted_diams_400msMUA.npy']);
writeNPY(stimParams.TF, [sDir,'TF_400msMUA.npy']);
writeNPY(stimParams.SF, [sDir,'SF_400msMUA.npy']);
writeNPY(contrast, [sDir,'contrast_400msMUA.npy']);

if plot_all
    for i = 1:size(spikeCnts_NoL,1)
        figure;
        errorbar(stimParams.VRF_Diam, mean(squeeze(spikeCnts_NoL(i,:,:)),2), std(squeeze(spikeCnts_NoL(i,:,:)),[],2)./sqrt(size(spikeCnts_NoL,3)),'ko-')
        set(gca,'XScale','log')
    end
end

