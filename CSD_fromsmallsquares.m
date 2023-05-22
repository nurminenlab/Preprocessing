% this script takes takes in an unfiltered LFP signal and outputs a
% CSD profile

base_dir    = '/opt3/';
animal      = 'MM385/';

l_coff = 3;
h_coff = 100;
butt_order = 2;

R = 0.2;
h = 0.25;
sigma = 0.4;

load('arrayLayoutV_Probe.mat')

addpath(genpath('/home/lauri/matlab_functions/NPMK'));
addpath(genpath('/home/lauri/code/npy-matlab'));
addpath(genpath('/home/lauri/code/spikes'));
addpath('/home/lauri/matlab_functions/statisticals');
addpath('/home/lauri/matlab_functions/kCSD');

sDir = uigetdir([base_dir, animal]);
sDir = [sDir, '/'];
NEVf = dir([sDir,'*-001.mat']);
NEV = openNEV([sDir,NEVf.name])
%load([sDir,'*-001.mat']);
%load([sDir,NEVf.name])

fs  = NEV.MetaTags.TimeRes;

stparf = dir([sDir, 'stimParams*.mat']);
load([stparf(1).folder,'/',stparf(1).name])

% get pointer to dat-file
datf = dir(fullfile(sDir,'*.ns5.dat'));

% open continuous TTL channel record
a(1) = dir([sDir,'20*.ns4']);
DataStruct  = openNSxNew([a(1).folder,'/',a(1).name],'read','c:1');
Fs          = DataStruct.MetaTags.SamplingFreq;

% collect stimulus info
anal_duration = stimParams.on_time_ms/1000 + 0.1;
base_duration = 0.1;
stimOnCount    = round((stimParams.on_time_ms/1000)*Fs);
stimOffCount   = round((stimParams.off_time_ms/1000)*Fs);
stimTotalCount = stimOnCount+stimOffCount;
stimMatrix = stimParams.stimOrder;
nStims    = size(stimParams.stimOrder,1);
nTrials   = size(stimParams.stimOrder,2);
sampFreqNEV = double(fs);

% generate pulse sample times from continuous trigger channel
pulseTimes = PulseCounter(DataStruct.Data,20000);
stimTimes  = reshape(pulseTimes,2,[]);
stimTimes  = stimTimes(1,:)+stimOffCount;
% note that these times refer to the beginning of a block
stimTimesNEV   = round(stimTimes*(sampFreqNEV/Fs));

Nchan = 32;
NSampRead = round((anal_duration + 2*base_duration)*30e3);

DD = dir([sDir,'*.ns5']);
NS = openNSx([sDir, DD.name]);
Data = double(NS.Data(arrayLayout(1:24),:));

N = 10; %order
F60 = 60; %notch frequency to filter line noise
F120 = 120; %notch frequency to filter harmonics of the line noise
Q = 10; %quality factor
Hd60 = design(fdesign.notch('N,F0,Q',N,F60,Q,fs));
Hd120 = design(fdesign.notch('N,F0,Q',N,F120,Q,fs));
d = fdesign.bandpass('N,F3dB1,F3dB2',4,3,300,fs);
HdPB = design(d,'butter');

for k=1:size(Data,1)
    fprintf('Filtering electrode %d \n',k)
    tmp = filtfilt(Hd60.sosMatrix,Hd60.ScaleValues,Data(k,:));
    tmp = filtfilt(Hd120.sosMatrix,Hd120.ScaleValues,tmp);
    tmp = filtfilt(HdPB.sosMatrix,HdPB.ScaleValues,tmp);
    Data(k,:) = tmp;
end

LFP_store = zeros(size(stimMatrix,1),24,round(anal_duration*30e3+0.1*30e3),size(stimMatrix,2));

% construct butterworth filter
stim_nums = unique(stimMatrix);
for f = 1:length(stim_nums)
    for g = 1:size(stimMatrix,2)
        % find block position
        b_pos = find(stimMatrix(:,g) == f);
        % to avoid negative sample at first stimulus
        epoch_start = stimTimesNEV(g) + (b_pos-1)*30e3;
        LFP_store(f,:,:,g) = Data(:,round(epoch_start - 0.1*30e3):round(epoch_start + anal_duration*30e3)-1);
    end
end

LFP = mean(LFP_store,4);
elPos = 0:0.1:2.3;
X = [-0.1:0.01:2.4];
SigEnd = 21000;
OnsetLFP = 0.1*30e3;
scale = max(max(max(LFP)));

CSD      = zeros(size(LFP,1),length(X),SigEnd);
CSD_nrmd = zeros(size(LFP,1),length(X),SigEnd);
% calc CSD
for i = 1:size(LFP,1)
    k = kCSD1d(elPos, squeeze(LFP(i,:,:)), 'X', X, 'R', R, 'h', h, 'sigma', sigma);
    k.chooseLambda();
    k.estimate;
    CSD(i,:,:) = k.csdEst;
end

% normd CSD
for i = 1:size(CSD,1)
    for j = 1:size(CSD,2)
        fak = std(squeeze(CSD(i,j,1:OnsetLFP)));
        CSD_nrmd(i,j,:) = (CSD(i,j,:) - mean(squeeze(CSD(i,j,1:OnsetLFP))))./fak;
    end
end

maxCSD      = nanmax(abs([nanmax(nanmax(nanmax(CSD))), nanmin(nanmin(nanmin(CSD)))]));
maxCSD_nrmd = nanmax(abs([nanmax(nanmax(nanmax(CSD_nrmd))), nanmin(nanmin(nanmin(CSD_nrmd)))]));

for i = 1:size(LFP,1)
    LFP2 = squeeze(LFP(i,:,:));
    figure(1)

    % LFP traces
    h1 = subplot(1,3,1);
    plot([1:size(LFP2,2)],LFP2 - repmat(((0:size(LFP2,1)-1)*0.75*scale)',1,size(LFP2,2)))
    Ylim2=get(h1,'YLim');
    Xlim2=get(h1,'XLim');
    set(h1,'xtick',[OnsetLFP,OnsetLFP+1500]);
    set(h1,'xtickLabel',[0,50]);
    set(h1,'YTick',fliplr(-((-1:size(LFP2,1)-1)*0.75*scale)));
    Yax_tick_loc = get(h1,'YTick');
    set(h1,'YTicklabel',fliplr([-0.1, elPos]),'FontSize',6);
    hold on
    plot([OnsetLFP OnsetLFP],[Ylim2(1) Ylim2(2)])
    axis([Xlim2(1) SigEnd -24*0.75*scale 1*0.75*scale ])
    ylabel('Electrode depth')
    xlabel('Peristimulus time (ms)')

    % CSD profile
    h2 = subplot(1,3,2);
    imagesc(squeeze(CSD(i,:,:)),[-maxCSD, maxCSD])
    hold on
    Ylim = get(h2,'YLim');
    plot([OnsetLFP OnsetLFP],[Ylim(1) Ylim(2)],'k-','LineWidth',2);
    set(h2,'YTick',1:10:250)
    set(h2,'YTickLabel',[-0.1:0.1:2.4])
    set(h2,'xtick',[OnsetLFP,OnsetLFP+1500]);
    set(h2,'xtickLabel',[0,50]);
    colormap(jet)
    
    % normalized CSD
    h3 = subplot(1,3,3);
    inds = isnan(squeeze(CSD_nrmd(i,:,:)));
    CSD_nrmd(i,squeeze(inds)) = 0;
    if maxCSD_nrmd == 0
        imagesc(squeeze(CSD_nrmd(i,:,:)))
    else
        imagesc(squeeze(CSD_nrmd(i,:,:)),[-maxCSD_nrmd, maxCSD_nrmd])
    end

    hold on
    Ylim = get(h3,'YLim');
    plot([OnsetLFP OnsetLFP],[Ylim(1) Ylim(2)],'k-','LineWidth',2);
    set(h3,'YTick',1:10:250)
    set(h3,'YTickLabel',[-0.1:0.1:2.4])
    set(h3,'xtick',[OnsetLFP,OnsetLFP+1500]);
    set(h3,'xtickLabel',[0,50]);
    print(['CSD', num2str(i)], '-depsc')
    clf;
end

% make the CSD plot
for i = 1:size(LFP,1)
    LFP2 = squeeze(LFP(i,:,:));
    figure(2)
    % CSD profile
    h2 = subplot(sqrt(size(LFP,1)),sqrt(size(LFP,1)),i);
    imagesc(squeeze(CSD(i,:,:)),[-maxCSD, maxCSD])
    hold on
    Ylim = get(h2,'YLim');
    plot([OnsetLFP OnsetLFP],[Ylim(1) Ylim(2)],'k-','LineWidth',2);
    set(h2,'YTick',1:10:250)
    set(h2,'YTickLabel',[-0.1:0.1:2.4])
    set(h2,'xtick',[OnsetLFP,OnsetLFP+1500]);
    set(h2,'xtickLabel',[0,50]);
    colormap(jet)
    axis off
    set(h2,'XTick',[OnsetLFP,OnsetLFP+1500]);
    set(h2,'xtickLabel',[0,50]);
end

