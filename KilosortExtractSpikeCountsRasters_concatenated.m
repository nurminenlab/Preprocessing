function KilosortExtractSpikeCountsRasters_concatenated()
addpath(genpath('/home/lauri/matlab_functions/NPMK'));
addpath(genpath('/home/lauri/code/npy-matlab'));
addpath(genpath('/home/lauri/code/spikes'));
addpath('/home/lauri/matlab_functions/statisticals');

%close all
errb = 1;
clc;

base_dir    = '/opt3/';
animal      = 'MM390/';
penetration = 'P9/';

anal_offset   = 0.05;
anal_duration = 0.4;
base_duration = 0.4;
base_onset    = 0.1;
base_offset   = base_onset + anal_duration;

% choose folder where the concatenated files live
sDir = uigetdir('/opt3/');
sDir = [sDir, '/'];
cat_info = load(fullfile(sDir,'concatDataInfo.mat'));

%% spikes 
ksf = dir(fullfile(sDir,'rez.mat'));
kilorez = load([sDir,ksf(1).name]);
NEVfreq = kilorez.rez.ops.fs;
% load spike times 
sptf = dir(fullfile(sDir,'spike_times.npy'));
spt  = readNPY([sDir,sptf(1).name]);
% load cluster identities 
cluf = dir(fullfile(sDir,'spike_clusters.npy'));
clu  = readNPY([sDir,cluf(1).name]);

for cat_f = 1:length(cat_info.fileNames)

    stparf = dir([base_dir,animal,penetration,cat_info.fileNames{cat_f}(1:end-8),'/','stimParams*.mat']);
    load([stparf(1).folder,'/',stparf(1).name])

    % get pointer to dat-file
    datf = dir(fullfile(sDir,'*combined.ns5.dat'));

    % in ms defines how much before and after stimulus onset and offset
    % is taken for spike-rasters
    raster_trail = 100;
    raster_befre = 400;
    
    % open continuous TTL channel record
    a(1) = dir([base_dir,animal,penetration,cat_info.fileNames{cat_f}(1:end-8),'/','20*.ns4']);
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
    sampFreqNEV   = double(kilorez.rez.ops.fs);

    % because recordings were catenated
    if cat_f == 1;
        stimTimesNEV  = round(stimTimes*(sampFreqNEV/Fs));
    else
        stimTimesNEV  = cat_info.fileDurations(cat_f-1) +  round(stimTimes*(sampFreqNEV/Fs));
    end
    
    stimLengthNEV = round(stimOnCount*(sampFreqNEV/Fs));
    stimOffNEV    = round(stimOffCount*(sampFreqNEV/Fs));
    % Use time after last stimulus for blank response
    blnkTimes = stimTimesNEV(end,:) + stimLengthNEV + 0.1*30e3;

    % conf laser info
    Laser_Mat      = reshape(stimParams.laser.rnd_laser',2,(size(stimParams.laser.rnd_laser,2)*size(stimParams.laser.rnd_laser,1))/2);
    Laser_Mat(2,:) = [];

    N_units       = length(unique(clu));
    units         = unique(clu);
    
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
        if u == 1 & cat_f == 1
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

    % get spike info
    tempsUnW       = readNPY([sDir,'templates_unw.npy']);
    spikeTemplates = readNPY([sDir,'spike_templates.npy']);
    temps          = readNPY([sDir,'templates.npy']);
    WInv           = readNPY([sDir,'whitening_mat_inv.npy']);
    coords         = readNPY([sDir,'channel_positions.npy']);
    amplitudes     = readNPY([sDir,'amplitudes.npy']);

    % compute spike-widths in milliseconds using cortex-lab script
    [spikeWidths, tempWidths] = computeSpikeWidths(tempsUnW, spikeTemplates);
    mn_spikeWidths = NaN*ones(size(units));
    for u = 1:length(units)
        mn_spikeWidths(u) = abs(mean(spikeWidths(clu == units(u))./NEVfreq*1000));
    end

    % get waveforms and depths
    mn_spikeDepths = NaN*ones(size(units));
    [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = templatePositionsAmplitudes(temps, WInv, coords(:,2), spikeTemplates, amplitudes);
    for u = 1:length(units)
        mn_spikeDepths(u) = mean(spikeDepths(clu == units(u)));
    end

    % define some parameters for SNR calculation
    gwfparams.dataDir = sDir;
    gwfparams.dataType = 'int16';   % Data type of .dat file (this should be BP filtered)
    gwfparams.nCh      = 32;        % Number of channels that were streamed to disk in .dat file
    gwfparams.wfWin    = [-40 41];  % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf      = 10e3;      % Max number of waveforms per unit to pull out
    gwfparams.fshigh   = 300;
    gwfparams.fs       = 30e3;

    fileName = fullfile(sDir,datf.name);
    filenamestruct = dir(fileName);
    dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8'));
    gwfparams.nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes); 

    fprintf('High-pass filtering\n')
    filtred_data = batch_filter_trace(fileName,kilorez.rez.ops,'high',3,300);

    fprintf('Extracting SNR\n')
    SNR       = NaN*ones(size(units));
    waveForms = NaN*ones(size(units,1),10e3,82);
    mxmp_cnt  = NaN*ones(size(units));
    for u = 1:length(units)
        gwfparams.spikeTimes = spt(clu == units(u));
        gwfparams.spikeClusters = repmat(units(u), size(gwfparams.spikeTimes));
        [SNR(u),waveForms(u,:,:),mxmp_cnt(u)] = getSNR(gwfparams,filtred_data);
    end

    diams = unique(stimParams.Value)';

    if cat_f == 1;
        spikeCnts_L   = spikeCnts_L_tmp;
        spikeCnts_NoL = spikeCnts_NoL_tmp;
        baseLine      = baseLine_tmp;
        spkraster_L   = spkraster_L_tmp;
        spkraster_NoL = spkraster_NoL_tmp;
        contrast      = stimParams.Contrast*ones(size(spikeCnts_L_tmp,3),1);
    else
        spikeCnts_L   = cat(3,spikeCnts_L,spikeCnts_L_tmp);
        spikeCnts_NoL = cat(3,spikeCnts_NoL,spikeCnts_NoL_tmp);
        baseLine      = cat(3,baseLine,baseLine_tmp);
        spkraster_L   = cat(3,spkraster_L,spkraster_L_tmp);
        spkraster_NoL = cat(3,spkraster_NoL,spkraster_NoL_tmp);
        contrast      = cat(1,contrast,stimParams.Contrast*ones(size(spikeCnts_L_tmp,3),1));
    end

end

writeNPY(spikeCnts_L, [sDir,'Kilosorted_spkCnts_L_concat400ms.npy']);
writeNPY(spikeCnts_NoL, [sDir,'Kilosorted_spkCnts_NoL_concat400ms.npy']);
writeNPY(baseLine, [sDir,'Kilosorted_baseLine_concat400ms.npy']);
writeNPY(spkraster_L, [sDir,'Kilosorted_spkraster_L_concat400ms.npy']);
writeNPY(spkraster_NoL, [sDir,'Kilosorted_spkraster_NoL_concat400ms.npy']);
writeNPY(fit_params_L, [sDir,'Kilosorted_fitParams_L_concat400ms.npy']);
writeNPY(fit_params, [sDir,'Kilosorted_fitParams_NoL_concat400ms.npy']);
writeNPY(field_params_L, [sDir,'Kilosorted_fieldParams_L_concat400ms.npy']);
writeNPY(field_params, [sDir,'Kilosorted_fieldParams_NoL_concat400ms.npy']);
writeNPY(mn_spikeWidths, [sDir,'Kilosorted_spikeWidths_concat400ms.npy']);
writeNPY(mxmp_cnt./10, [sDir,'Kilosorted_spikeDepths_concat400ms.npy']);
writeNPY(SNR, [sDir,'Kilosorted_SNR_concat400ms.npy']);
writeNPY(waveForms, [sDir,'Kilosorted_waveForms_concat400ms.npy']);
writeNPY(mxmp_cnt./10,[sDir,'Kilosorted_maxamp_contact_concat400ms.npy']);
writeNPY(diams, [sDir,'Kilosorted_diams_concat400ms.npy']);
writeNPY(stimParams.TF, [sDir,'TF_concat400ms.npy']);
writeNPY(stimParams.SF, [sDir,'SF_concat400ms.npy']);
writeNPY(contrast, [sDir,'contrast_concat400ms.npy']);
