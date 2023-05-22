function ExtractSpikeCountsRasters()
addpath(genpath('/home/lauri/matlab_functions/NPMK'));
addpath(genpath('/home/lauri/code/npy-matlab'));
addpath(genpath('/home/lauri/code/spikes'));
addpath('/home/lauri/matlab_functions/statisticals');

%close all
errb = 1;
clc;

base_dir    = '/opt3/';
animal      = 'MM378/';
penetration = 'P2/';

anal_offset   = 0.05;
anal_duration = 0.4;
base_duration = 0.4;
base_onset    = 0.1;
base_offset   = base_onset + anal_duration;

% choose folder where the concatenated files live
sDir = uigetdir('/opt3/');
sDir = [sDir, '/'];

%% spikes 
NEV = openNEV([sDir,'thresholdedNEV.nev-rthr.nev'])
fs  = NEV.MetaTags.TimeRes;
spt = NEV.Data.Spikes.TimeStamp;
clu = NEV.Data.Spikes.Electrode;

stparf = dir([sDir,'stimParams*.mat']);
load([stparf(1).folder,'/',stparf(1).name])

% get pointer to dat-file
datf = dir(fullfile(sDir,'*-001.ns5.dat'));

% in ms defines how much before and after stimulus onset and offset
% is taken for spike-rasters
raster_trail = 100;
raster_befre = 400;
    
% open continuous TTL channel record
a(1) = dir([sDir,'/','20*.ns4']);
DataStruct  = openNSxNew([a(1).folder,'/',a(1).name],'read','c:1');
Fs          = DataStruct.MetaTags.SamplingFreq;
% generate pulse sample times from continuous trigger channel
pulseTimes  = PulseCounter(DataStruct.Data,20000);

%% collect stimulus info
stimOnCount    = round((stimParams.on_time_ms/1000)*Fs);
stimOffCount   = round((stimParams.off_time_ms/1000)*Fs);
stimTotalCount = stimOnCount+stimOffCount;
stimMatrix_tmp = stimParams.stimOrder;
stimMatrix     = ones(size(stimMatrix_tmp,1)/2,size(stimMatrix_tmp,2))*NaN;
% every stimulus is run twice in succession with and without the laser
for i = 1:size(stimMatrix_tmp,2)
    stimMatrix(:,i)     = unique(stimMatrix_tmp(:,i),'stable');
end
nStims    = size(stimParams.stimOrder,1);
nTrials   = size(stimParams.stimOrder,2);
% div by 2 because just every second stimulus produces a trigger
stimTimes = reshape(pulseTimes,nStims/2,nTrials); 

%% calculate number of spikes
sampFreqNEV   = double(fs);
stimTimesNEV  = round(stimTimes*(sampFreqNEV/Fs));
    
stimLengthNEV = round(stimOnCount*(sampFreqNEV/Fs));
stimOffNEV    = round(stimOffCount*(sampFreqNEV/Fs));
% Use time after last stimulus for blank response
blnkTimes = stimTimesNEV(end,:) + stimLengthNEV + 0.1*30e3;

% conf laser info
Laser_Mat      = reshape(stimParams.laser.rnd_laser',2,(size(stimParams.laser.rnd_laser,2)*size(stimParams.laser.rnd_laser,1))/2);
Laser_Mat(2,:) = [];

N_units       = length(unique(NEV.Data.Spikes.Electrode));
units         = unique(NEV.Data.Spikes.Electrode);
    
spikeCnts_L_tmp   = NaN*ones(N_units, nStims/2, nTrials);
spikeCnts_NoL_tmp = NaN*ones(N_units, nStims/2, nTrials);
baseLine_tmp      = NaN*ones(N_units, nStims/2, nTrials);
spkraster_L_tmp   = NaN*ones(N_units, nStims, nTrials, length(-raster_befre:1:stimParams.on_time_ms+raster_trail)-1);
spkraster_NoL_tmp = NaN*ones(N_units, nStims, nTrials, length(-raster_befre:1:stimParams.on_time_ms+raster_trail)-1);

for u = 1:N_units
    spikeTimeMat = double(spt(units(u) == clu));
    
    % collate spike-counts depending on weather the laser was on or off
    for f = 1:nStims/2
        idx = find(stimMatrix==f);        
        for g = 1:length(idx)

            % laser was on during the first stimulus of the pair
            if Laser_Mat(idx(g)) == 1
                % raster
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) - base_duration*30e3) & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV + base_duration*30e3);
                spk_tmp = (double(spikeTimeMat(spk_tmp) - stimTimesNEV(idx(g))))./30e3*1000;
                spkraster_L_tmp(u,f,g,:) = histcounts(spk_tmp,[-raster_befre:1:stimParams.on_time_ms+raster_trail]);
                % count
                spikeCnts_L_tmp(u,f,g)   = sum(spikeTimeMat>stimTimesNEV(idx(g))+anal_offset*30e3 & spikeTimeMat<stimTimesNEV(idx(g))+anal_offset*30e3+anal_duration*30e3);

                % DON
                % raster
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV - base_duration*30e3) & (spikeTimeMat<stimTimesNEV(idx(g))+ 2*stimLengthNEV + stimOffNEV + base_duration*30e3));
                spk_tmp = (double(spikeTimeMat(spk_tmp) - (stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV)))./30e3*1000;
                spkraster_NoL_tmp(u,f,g,:) = histcounts(spk_tmp,[-raster_befre:1:stimParams.on_time_ms+raster_trail]);
                % count
                spikeCnts_NoL_tmp(u,f,g) = sum(spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV+anal_offset*30e3) & spikeTimeMat<stimTimesNEV(idx(g))+ anal_offset*30e3+anal_duration*30e3+stimLengthNEV + stimOffNEV);

            % laser was on during the second stimulus of the pair
            else
                % raster
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) - base_duration*30e3) & (spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV + base_duration*30e3));
                spk_tmp = (double(spikeTimeMat(spk_tmp) - stimTimesNEV(idx(g))))./30e3*1000;
                spkraster_NoL_tmp(u,f,g,:) = histcounts(spk_tmp,[-raster_befre:1:stimParams.on_time_ms+raster_trail]);
                % count
                spikeCnts_NoL_tmp(u,f,g)   = sum(spikeTimeMat>stimTimesNEV(idx(g))+anal_offset*30e3 & spikeTimeMat<stimTimesNEV(idx(g))+anal_offset*30e3+anal_duration*30e3);

                % raster
                spk_tmp = (spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV - base_duration*30e3) & spikeTimeMat<stimTimesNEV(idx(g))+ 2*stimLengthNEV + stimOffNEV + base_duration*30e3);
                spk_tmp = (double(spikeTimeMat(spk_tmp) - (stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV)))./30e3*1000;
                spkraster_L_tmp(u,f,g,:) = histcounts(spk_tmp,[-raster_befre:1:stimParams.on_time_ms+raster_trail]);
                % count
                spikeCnts_L_tmp(u,f,g)   = sum(spikeTimeMat>(stimTimesNEV(idx(g)) + stimLengthNEV + stimOffNEV+anal_offset*30e3) & spikeTimeMat<stimTimesNEV(idx(g))+anal_offset*30e3+anal_duration*30e3+stimLengthNEV + stimOffNEV);

            end        
            baseLine_tmp(u,f,g) = sum(spikeTimeMat > stimTimesNEV(idx(g))+stimLengthNEV+base_onset*30e3 & spikeTimeMat<stimTimesNEV(idx(g))+stimLengthNEV+base_offset*30e3);
        end
    end

    fprintf('Now fitting unit %d\n', u)
    if u == 1
        f = fit_sizetuning(mean(spikeCnts_NoL_tmp(u,:,:),3),unique(stimParams.Value)');
        fit_params     = f.final_params;
        field_params   = f.field_params;

        f = fit_sizetuning(mean(spikeCnts_L_tmp(u,:,:),3),unique(stimParams.Value)');
        fit_params_L   = f.final_params;
        field_params_L = f.field_params;
    else
        f = fit_sizetuning(mean(spikeCnts_NoL_tmp(u,:,:),3),unique(stimParams.Value)');
        fit_params     = [fit_params; f.final_params];
        field_params   = [field_params; f.field_params];
        
        f = fit_sizetuning(mean(spikeCnts_L_tmp(u,:,:),3),unique(stimParams.Value)');
        fit_params_L = [fit_params_L; f.final_params];
        field_params_L = [field_params_L; f.field_params];
    end

end

mn_spikeWidths = NaN*ones(size(units));

SNR       = NaN*ones(size(units));
waveForms = NaN*ones(size(units,1),10e3,82);
mxmp_cnt  = NaN*ones(size(units));
mn_spikeDepths = NaN*ones(size(units));
    
diams = unique(stimParams.Value)';

spikeCnts_L   = spikeCnts_L_tmp;
spikeCnts_NoL = spikeCnts_NoL_tmp;
baseLine      = baseLine_tmp;
spkraster_L   = spkraster_L_tmp;
spkraster_NoL = spkraster_NoL_tmp;
contrast      = stimParams.Contrast*ones(size(spikeCnts_L_tmp,3),1);

writeNPY(spikeCnts_L, [sDir,'Kilosorted_spkCnts_L_concat_MUA400.npy']);
writeNPY(spikeCnts_NoL, [sDir,'Kilosorted_spkCnts_NoL_concat_MUA400.npy']);
writeNPY(baseLine, [sDir,'Kilosorted_baseLine_concat_MUA400.npy']);
writeNPY(spkraster_L, [sDir,'Kilosorted_spkraster_L_concat_MUA400.npy']);
writeNPY(spkraster_NoL, [sDir,'Kilosorted_spkraster_NoL_concat_MUA400.npy']);
writeNPY(fit_params_L, [sDir,'Kilosorted_fitParams_L_concat_MUA400.npy']);
writeNPY(fit_params, [sDir,'Kilosorted_fitParams_NoL_concat_MUA400.npy']);
writeNPY(field_params_L, [sDir,'Kilosorted_fieldParams_L_concat_MUA400.npy']);
writeNPY(field_params, [sDir,'Kilosorted_fieldParams_NoL_concat_MUA400.npy']);
writeNPY(mn_spikeWidths, [sDir,'Kilosorted_spikeWidths_concat_MUA400.npy']);
% valid depth because batch_filter_trace.m reorganizes electrodes
writeNPY((1:size(spikeCnts_NoL,1))/10, [sDir,'Kilosorted_spikeDepths_concat_MUA400.npy']);
writeNPY(SNR, [sDir,'Kilosorted_SNR_concat_MUA400.npy']);
writeNPY(waveForms, [sDir,'Kilosorted_waveForms_concat_MUA400.npy']);
writeNPY(mxmp_cnt,[sDir,'Kilosorted_maxamp_contact_concat_MUA400.npy']);
writeNPY(diams, [sDir,'Kilosorted_diams_concat_MUA400.npy']);
writeNPY(stimParams.TF, [sDir,'TF_concat_MUA400.npy']);
writeNPY(stimParams.SF, [sDir,'SF_concat_MUA400.npy']);
writeNPY(contrast, [sDir,'contrast_concat_MUA400.npy']);

