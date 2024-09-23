% the aim of this script is to prepare data for the Master Class
% Demonstration in Spain

clear all; close all; clc

DataPFC = load('PFC_100_CH63_0.continuous.mat');
DataHPC = load('HPC_100_CH8_0.continuous.mat');
DataPFC = DataPFC.PFC;
DataHPC = DataHPC.HPC;

SampF = 1000; 
NDownLFP = 2; % how many times you want to reduce the sampling of the data 
FsDownLFP = SampF/NDownLFP;
lfpPFCDown = decimate(DataPFC,NDownLFP,'FIR');
lfpHPCDown = decimate(DataHPC,NDownLFP,'FIR');
timVect1 = linspace(0,numel(lfpPFCDown)/FsDownLFP,numel(lfpPFCDown));
%%
% plot the results to inspect them
figure
plot(timVect1, lfpHPCDown+3000,'r')
hold on 
plot(timVect1, lfpPFCDown, 'b')
legend('HPC', 'PFC')
xlabel('Time [s]')
ylabel('Amplitude [\mu V]')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
%print(gcf,'-dpdf', 'signal.pdf')
%% Zoom_in specific parts of the signal - theta
figure
plot(timVect1, lfpHPCDown,'r')
xlim([2300 2305])
legend('HPC', 'PFC')
xlabel('Time [s]')
ylabel('Amplitude [\mu V]')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
%print(gcf,'-dpdf', 'theta.pdf')

%% Zoom-In Delta Oscillations
figure
plot(timVect1, lfpPFCDown,'b')
xlim([2040 2045])
legend('HPC', 'PFC')
xlabel('Time [s]')
ylabel('Amplitude [\mu V]')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
%print(gcf,'-dpdf', 'delta.pdf')
%%
downSampleFreq = 2;
ampThresh = [4 8];
timeWinThresh = [1 0.1];
[lfpHPCDownN,artefactIndsHPC] = removeArtefacts(lfpHPCDown,downSampleFreq,ampThresh, ...
    timeWinThresh);
[lfpPFCDownN,artefactIndsPFC] = removeArtefacts(lfpPFCDown,downSampleFreq,ampThresh, ...
    timeWinThresh);
figure
plot(timVect1, lfpHPCDownN+2000,'r')
hold on 
plot(timVect1, lfpPFCDownN, 'b')
legend('HPC', 'PFC')
xlabel('Time [s]')
ylabel('Amplitude [\mu V]')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
%print(gcf,'-dpdf', 'signal_artefact_free.pdf')
disp('Finished')
%%
% 1. the parameters for the computation 
samplingFrequency = 5;
Fs = SampF; 
binScootS = 1 ./ samplingFrequency;
binScootSamps = round(Fs*binScootS);

% 2. You filter the signal in the high frequency range
% Filter first in high frequency band to remove low-freq physiologically
% correlated LFPs (e.g., theta, delta, SPWs, etc.)

maxfreqband = floor(max([500 Fs/2]));
xcorr_freqband = [275 300 maxfreqband-25 maxfreqband]; % Hz
filteredPFC  = filtsig_in(lfpPFCDownN, Fs, xcorr_freqband);
filteredHPC  = filtsig_in(lfpHPCDownN, Fs, xcorr_freqband);

% 3. The relevant parameters important for further correlation analysis 
xcorr_window_samps = round(binScootS*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples
timestamps = (1+xcorr_window_inds(end)):binScootSamps:(size(filteredPFC,1)-xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
counter = 1;
c1 = [];
c2 = [];
binind = 0;
binindstart = 1;

for i = 1: numel(timestamps)
    binind = binind+1;
    s1 =filteredPFC(timestamps(i) + xcorr_window_inds);
    s2 =filteredHPC(timestamps(i)+ xcorr_window_inds);
    c1 = cat(2,c1,s1);
    c2 = cat(2,c2,s2);
    binindend = binind;
    tmp = corr(c1,c2);
    tmp = diag(tmp);
    EMGCorr(binindstart:binindend) = EMGCorr(binindstart:binindend) + tmp;
    c1 = [];
    c2 = [];
    binindstart = binind+1;
end
EMGCorr = EMGCorr/(2*(2-1)/2);
EMGNorm = bz_NormToRange(EMGCorr,[0 1]);

% 4. Making the final structure for the EMG data 
EMGFromLFP.timestamps = timestamps'./Fs;
EMGFromLFP.data = EMGCorr;
EMGFromLFP.Norm = EMGNorm;
EMGFromLFP.channels = 'HPC and PFC';
EMGFromLFP.detectorName = 'bz_EMGFromLFP';
EMGFromLFP.samplingFrequency = samplingFrequency; 

% 5. Smoothing the outcome of the EMG 
smoothfact = 20;
dtEMG = 1/EMGFromLFP.samplingFrequency;
EMGFromLFP.smoothed = smooth(EMGFromLFP.data,smoothfact/dtEMG,'moving');
disp('Finished')
%%
figure 
plot(EMGFromLFP.timestamps,EMGFromLFP.smoothed ...
    ,'LineWidth',2)
xlabel('Time [s]')
ylabel('A. U.')
title('EMG-like signal 5 Hz SF')
legend boxoff
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
box off  
set(gcf,'Color','w')
%print(gcf,'-dpdf', 'emg_like_signal.pdf')

%%
% 1. The important parameters for the analysis 
window = 2;
numhistbins = 21;
histbins = linspace(0,1,numhistbins);
numfreqs = 100;
swFFTfreqs = logspace(0,2,numfreqs);
noverlap = window-1; %Updated to speed up (don't need to sample at fine time resolution for channel selection)
SWweights = 'PSS';
SWWeightsName = 'PSS';
SWfreqlist = logspace(0.5,2,numfreqs); %should get this from bz_PowerSpectrumSlope...

%For Theta Calculation
f_all = [1 20];
f_theta = [5 10];
thFFTfreqs = logspace(log10(f_all(1)),log10(f_all(2)),numfreqs);
numSWChannels = 1;
numThetaChannels = 1;
swhists = zeros(numhistbins,numSWChannels);
dipSW = zeros(numSWChannels,1);
THhist = zeros(numhistbins,numThetaChannels);
THmeanspec = zeros(numfreqs,numThetaChannels);
peakTH = zeros(numThetaChannels,1);
IRASA = 'TRUE';
ThIRASA = 'TRUE';
% ===================================================== %
% 2. Calculating the slow wave metric using the slope of the power spectrum

display('Calculating SW Mertric using Power Spectrum Slope')
%Put the LFP in the right structure format
lfp.data = lfpPFCDownN;
lfp.timestamps = timVect1;
lfp.samplingRate = FsDownLFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,'frange',[1 90],'IRASA',IRASA);
broadbandSlowWave = -specslope.data; %So NREM is higher as opposed to lower
t_clus = specslope.timestamps;
swFFTfreqs = specslope.freqs';
specdt = 1./specslope.samplingRate;
swFFTspec = 10.^spec.amp'; %To reverse log10 in bz_PowerSpectrumSlope
IRASAsmooth = spec.IRASAsmooth';
IRASAintercept = specslope.intercept;
IRASAslope = specslope.data;
    
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(log10(swFFTspec)','modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes = find(totz>3);

%Smooth and 0-1 normalize
broadbandSlowWave(badtimes) = nan;
broadbandSlowWave = smooth(broadbandSlowWave,smoothfact./specdt);
%% 5. Calculate theta power 
display('Calculating Theta Metric above PSS')
%Put the LFP in the right structure format
lfp.data = lfpHPCDownN;
lfp.timestamps = timVect1;
lfp.samplingRate = FsDownLFP;
%Calculate PSS
[specslope,spec] = bz_PowerSpectrumSlope(lfp,window,window-noverlap,...
        'nfreqs',200,'frange',f_all,'IRASA',ThIRASA);
t_thclu = specslope.timestamps;
specdt = 1./specslope.samplingRate;
thFFTspec = specslope.resid';
thFFTspec(thFFTspec<0)=0;
IRASAsmooth_th = spec.IRASAsmooth';
thFFTspec_raw = 10.^spec.amp';
    
% Remove transients before calculating SW histogram
zFFTspec = NormToInt(spec.amp,'modZ');
totz = NormToInt(abs(sum(zFFTspec,2)),'modZ');
badtimes_TH = find(totz>3);
    
thFFTfreqs = specslope.freqs';
thfreqs = (thFFTfreqs>=f_theta(1) & thFFTfreqs<=f_theta(2));
thratio = max((thFFTspec(thfreqs,:)),[],1);
thratio(badtimes_TH) = nan;     %Remove transients
thratio = smooth(thratio,smoothfact./specdt);
 % normalize your ration to the interval between 0 and 1 
broadbandSlowWave = bz_NormToRange(broadbandSlowWave,[0 1]);
thratio = bz_NormToRange(thratio,[0 1]);

%% 7. Making a plot that includes both ratios 
figure 
plot(t_clus,broadbandSlowWave+3,'b','LineWidth',2)
hold on
plot(t_clus,thratio,'r','LineWidth',2)
legend('Slow wave','Theta ratio','Location','best')
xlabel('Time [s]')
ylabel('a. u.')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
box off 
set(gcf,'Color','w')
f=gcf;
set(f,'Units','centimeters', 'Color','w')
f.Position(3)=26;
f.Position(4)=13;
f.PaperSize=[25 12];
%print(gcf,'-dpdf', 'theta_delta_oscillations3.pdf')

%% 8. Working on the EMG signal 
EMGN = EMGFromLFP;
%remove any t_clus before/after t_emg
prEMGtime = t_clus<EMGN.timestamps(1) | t_clus>EMGN.timestamps(end);
broadbandSlowWave(prEMGtime) = []; thratio(prEMGtime) = []; t_clus(prEMGtime) = [];

%interpolate to FFT time points;
EMG_interpolate = interp1(EMGN.timestamps,EMGN.smoothed,t_clus,'nearest');
disp('Finished handling the EMG data')
%%
% Divide PC1 for SWS
%DL Note: should replace all of this with calls to bz_BimodalThresh
numpeaks = 1;
numbins = 12;
%numbins = 12; %for Poster...
while numpeaks ~=2
    [swhist,swhistbins]= hist(broadbandSlowWave,numbins);
    
    [PKS,LOCS] = findpeaks_SleepScore(swhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

SWdiptest = bz_hartigansdiptest(broadbandSlowWave);
betweenpeaks = swhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-swhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

swthresh = betweenpeaks(diploc);


 
%SWS time points
NREMtimes = (broadbandSlowWave >swthresh);

figure
hold on
bar(swhistbins(swhistbins>swthresh),swhist(swhistbins>swthresh),'FaceColor','b','barwidth',0.9,'linewidth',1)
bar(swhistbins(swhistbins<=swthresh),swhist(swhistbins<=swthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
xline(swthresh,'k--','LineWidth',2)
ylabel('Count')
xlabel('Normalized Amplitude')
title('Histogram Delta')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
box off 
set(gcf,'Color','w')
f=gcf;
set(f,'Units','centimeters', 'Color','w')
f.Position(3)=10.5;
f.Position(4)=10.5;
f.PaperSize=[10 10];
%print(gcf,'-dpdf', 'histogram_delta_oscillations3.pdf')
%% theta histogram
numpeaks = 1;
numbins = 12;
if sum(isnan(EMG_interpolate))>0
   error('EMG seems to have NaN values...') 
end

while numpeaks ~=2
    %[EMGhist,EMGhistbins]= hist(EMG(NREMtimes==0),numbins); 
    [EMGhist,EMGhistbins]= hist(EMG_interpolate,numbins); %changed back 6/18/19

    [PKS,LOCS] = findpeaks_SleepScore([0 EMGhist],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    
    if numpeaks ==100
        display('Something is wrong with your EMG')
        return
    end
end

EMGdiptest = bz_hartigansdiptest(EMG_interpolate);
betweenpeaks = EMGhistbins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks_SleepScore(-EMGhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

EMGthresh = betweenpeaks(diploc);

MOVtimes = (broadbandSlowWave(:)<swthresh & EMG_interpolate(:)>EMGthresh);

figure
hold on
bar(EMGhistbins(EMGhistbins>EMGthresh),EMGhist(EMGhistbins>EMGthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
bar(EMGhistbins(EMGhistbins<=EMGthresh),EMGhist(EMGhistbins<=EMGthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
xline(EMGthresh,'k--','LineWidth',2)
ylabel('Count')
xlabel('Normalized Amplitude')
title('Histogram EMG-Like Signal')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold','FontName','Times')
box off 
set(gcf,'Color','w')
f=gcf;
set(f,'Units','centimeters', 'Color','w')
f.Position(3)=10.5;
f.Position(4)=10.5;
f.PaperSize=[10 10];
%print(gcf,'-dpdf', 'histogram_emg_oscillations3.pdf')

%% 11. Then Divide Theta (During NonMoving)
numpeaks = 1;
numbins = 12;
while numpeaks ~=2 && numbins <=25
    %[THhist,THhistbins]= hist(thratio(SWStimes==0 & MOVtimes==0),numbins);
    [THhist,THhistbins]= hist(thratio(MOVtimes==0),numbins);

    [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
    LOCS = sort(LOCS);
    numbins = numbins+1;
    numpeaks = length(LOCS);
end

THdiptest = bz_hartigansdiptest(thratio(MOVtimes==0));

if numpeaks ~= 2
	display('No bimodal dip found in theta. Trying to exclude NREM...')

    numbins = 12;
    %numbins = 15; %for Poster...
    while numpeaks ~=2 && numbins <=25
        [THhist,THhistbins]= hist(thratio(NREMtimes==0 & MOVtimes==0),numbins);

        [PKS,LOCS] = findpeaks_SleepScore(THhist,'NPeaks',2,'SortStr','descend');
        LOCS = sort(LOCS);
        numbins = numbins+1;
        numpeaks = length(LOCS);

    end
    
	try
        THdiptest = bz_hartigansdiptest(thratio(NREMtimes==0 & MOVtimes==0));
    catch
	end
end

if length(PKS)==2
    betweenpeaks = THhistbins(LOCS(1):LOCS(2));
    [dip,diploc] = findpeaks_SleepScore(-THhist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

    THthresh = betweenpeaks(diploc);

    REMtimes = (broadbandSlowWave<swthresh & EMG_interpolate<EMGthresh & thratio>THthresh);
else
    display('No bimodal dip found in theta. Use TheStateEditor to manually select your threshold (hotkey: A)')
    THthresh = 0;
%     REMtimes =(broadbandSlowWave<swthresh & EMG<EMGthresh);
end

figure
hold on
bar(THhistbins(THhistbins>=THthresh),THhist(THhistbins>=THthresh),'FaceColor','r','barwidth',0.9,'linewidth',1)
bar(THhistbins(THhistbins<THthresh),THhist(THhistbins<THthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
xline(THthresh,'k--','LineWidth',2)
ylabel('Count')
xlabel('Normalized Amplitude')
title('Histogram Theta')
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
box off 
set(gcf,'Color','w')
f=gcf;
set(f,'Units','centimeters', 'Color','w')
f.Position(3)=10.5;
f.Position(4)=10.5;
f.PaperSize=[10 10];
%print(gcf,'-dpdf', 'histogram_theta_oscillations3.pdf')

%% 11. Now I have to create an output structure to include all these variables 
histsandthreshs = v2struct(swhist,swhistbins,swthresh,EMGhist,EMGhistbins,...
    EMGthresh,THhist,THhistbins,THthresh);
% Ouput Structure: StateScoreMetrics
% LFPparams = SleepScoreLFP.params;
WindowParams.window = window;
WindowParams.smoothwin = smoothfact;
% THchanID = SleepScoreLFP.THchanID; SWchanID = SleepScoreLFP.SWchanID;
% SleepScoreMetrics = v2struct(broadbandSlowWave,thratio,EMG,...
%     t_clus,badtimes,badtimes_TH,histsandthreshs,LFPparams,WindowParams,THchanID,SWchanID,...
%     recordingname,THdiptest,EMGdiptest,SWdiptest);
EMG = EMG_interpolate;
SleepScoreMetrics = v2struct(broadbandSlowWave,thratio,EMG,...
    t_clus,badtimes,badtimes_TH,histsandthreshs,WindowParams,...
    THdiptest,EMGdiptest,SWdiptest);
StatePlotMaterials = v2struct(t_clus,swFFTfreqs,swFFTspec,thFFTfreqs,thFFTspec);
if exist('IRASAsmooth_th','var')
    StatePlotMaterials.IRASAsmooth_th = IRASAsmooth_th;
    StatePlotMaterials.thFFTspec_raw = thFFTspec_raw;

end
if exist('IRASAsmooth','var')
    StatePlotMaterials.IRASAsmooth = IRASAsmooth;
    StatePlotMaterials.IRASAintercept = IRASAintercept;
    StatePlotMaterials.IRASAslope = IRASAslope;
end
%% 12 Determining the sleep states 
minSWSsecs = 6;
minWnexttoREMsecs = 6;
minWinREMsecs = 6;       
minREMinWsecs = 6;
minREMsecs = 6;
minWAKEsecs = 6;
MinTimeWindowParms = v2struct(minSWSsecs,minWnexttoREMsecs,minWinREMsecs,...
        minREMinWsecs,minREMsecs,minWAKEsecs);
[ints, idx, MinTimeWindowParms] = ClusterStates_DetermineStates(...
    SleepScoreMetrics,MinTimeWindowParms,histsandthreshs);
[zFFTspec,mu,sig] = zscore(log10(swFFTspec)');
if sum(isinf(log10(thFFTspec(:))))==0
    [~,mu_th,sig_th] = zscore(log10(thFFTspec)');
else %For Theta over PSS (ThIRASA)
    [~,mu_th,sig_th] = zscore((thFFTspec)');
end
viewwin  =[t_clus(1) t_clus(end)];
states = {ints.WAKEstate,ints.NREMstate,ints.REMstate};
figure 
subplot(8,1,[1:2])
imagesc(t_clus,log2(swFFTfreqs),log10(swFFTspec))
axis xy
set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
caxis([3.5 6.5])
caxis([min(mu)-2*max(sig) max(mu)+2*max(sig)])
xlim(viewwin)
colorbar('east')
ylim([log2(swFFTfreqs(1)) log2(swFFTfreqs(end))+0.2])
set(gca,'XTickLabel',{})
ylabel({'swLFP','f (Hz)'})
title('State Scoring Results');
subplot(8,1,3)
imagesc(t_clus,log2(thFFTfreqs),log10(thFFTspec))
axis xy
set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
%caxis([3.5 6.5])
caxis([min(mu_th)-2*max(sig_th) max(mu_th)+2*max(sig_th)])
xlim(viewwin)
%colorbar('east')
ylim([log2(thFFTfreqs(1)) log2(thFFTfreqs(end))+0.2])
ylabel({'thLFP','f (Hz)'})
set(gca,'XTickLabel',{})
subplot(8,1,4)
%plot(t_clus,-IDX,'LineWidth',2)
hold on
plot(states{1}',-1*ones(size(states{1}))','k','LineWidth',8)
plot(states{2}',-2*ones(size(states{2}))','b','LineWidth',8)
plot(states{3}',-3*ones(size(states{3}))','r','LineWidth',8)
xlim(viewwin)
ylim([-4 0])
set(gca,'YTick',[-3:-1])
set(gca,'YTickLabel',{'REM','SWS','Wake/MA'})
set(gca,'XTickLabel',{})
subplot(6,1,4)
hold on
plot(t_clus,broadbandSlowWave,'k','LineWidth',2)
%plot(synchtimes',thresh*ones(size(synchtimes))','r')
ylabel('SW')
box on
ylim([0 1])
xlim(viewwin)
set(gca,'XTickLabel',{})
subplot(6,1,5)
hold on
plot(t_clus,thratio,'k','LineWidth',2)
%plot(synchtimes',thresh*ones(size(synchtimes))','r')
ylabel('Theta')
box on
ylim([0 1])
xlim(viewwin)
set(gca,'XTickLabel',{})
subplot(6,1,6)
hold on
plot(t_clus,EMG,'k','LineWidth',2)
%plot(synchtimes',thresh*ones(size(synchtimes))','r')
ylabel('EMG')
box on
ylim([0 1])
xlim(viewwin)
xlabel('t (s)')
% set(gca,'FontSize',10,'LineWidth',1.5,'FontWeight','bold')
set(gcf,'Color','w')
% saveas(clusterfig,[figloc,recordingname,'_SSResults'],'jpeg')

f=gcf;
set(f,'Units','centimeters', 'Color','w')
f.Position(3)=30.5;
f.Position(4)=30.5;
f.PaperSize=[30 30];
%print(gcf,'-dpdf', 'sleep_state_scoring.pdf')
%%
%% 15 Split REM and arousal 
IDX_struct = bz_INTtoIDX(ints);
IDXNew = interp1(IDX_struct.timestamps,IDX_struct.states,t_clus,'nearest');
NREMtimes = (broadbandSlowWave >swthresh);
subplot(3,2,1)
hold on
bar(swhistbins(swhistbins>swthresh),swhist(swhistbins>swthresh),'FaceColor','b','barwidth',0.9,'linewidth',1)
bar(swhistbins(swhistbins<=swthresh),swhist(swhistbins<=swthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
plot([swthresh swthresh],[0 max(swhist)],'r','LineWidth',1)
xlabel('Broadband Slow Wave')
title('Step 1: Broadband for NREM')
subplot(3,2,3)
hold on
bar(EMGhistbins(EMGhistbins>EMGthresh),EMGhist(EMGhistbins>EMGthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
bar(EMGhistbins(EMGhistbins<=EMGthresh),EMGhist(EMGhistbins<=EMGthresh),'FaceColor',0.9*[1 1 1],'barwidth',0.9,'linewidth',1)
plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r','LineWidth',1)
xlabel('EMG')
title('Step 2: EMG for Muscle Tone')
subplot(3,2,5)
hold on
bar(THhistbins(THhistbins>=THthresh),THhist(THhistbins>=THthresh),'FaceColor','r','barwidth',0.9,'linewidth',1)
bar(THhistbins(THhistbins<THthresh),THhist(THhistbins<THthresh),'FaceColor','k','barwidth',0.9,'linewidth',1)
plot([THthresh THthresh],[0 max(THhist)],'r','LineWidth',1)
xlabel('Theta')
title('Step 3: Theta for REM')
subplot(2,2,2)
plot(broadbandSlowWave(IDXNew==2,1),EMG(IDXNew==2),'b.','markersize',5)
hold on
plot(broadbandSlowWave(EMG>EMGthresh & IDXNew==1,1),EMG(EMG>EMGthresh & IDXNew==1),'k.','markersize',5)
plot(broadbandSlowWave(EMG<EMGthresh & IDXNew==1|IDXNew==3,1),EMG(EMG<EMGthresh & IDXNew==1|IDXNew==3),...
            '.','Color',0.8*[1 1 1],'markersize',5)
plot(swthresh*[1 1],get(gca,'ylim'),'r','LineWidth',1)
plot(swthresh*[0 1],EMGthresh*[1 1],'r','LineWidth',1)
xlabel('Broadband SW');ylabel('EMG')
subplot(2,2,4)
plot(thratio(NREMtimes==0 & IDXNew==1,1),EMG(NREMtimes==0 & IDXNew==1,1),'k.','markersize',5)
hold on
plot(thratio(NREMtimes==0 & IDXNew==3,1),EMG(NREMtimes==0 & IDXNew==3,1),'r.','markersize',5)
xlabel('Narrowband Theta');ylabel('EMG')
plot(THthresh*[1 1],EMGthresh*[0 1],'r','LineWidth',1)
plot([0 1],EMGthresh*[1 1],'r','LineWidth',1)

f=gcf;
set(f,'Units','centimeters', 'Color','w')
f.Position(3)=20.5;
f.Position(4)=20.5;
f.PaperSize=[20 20];
%print(gcf,'-dpdf', 'allHistoTogether.pdf')
%%
%% 16 Figure: Clustering
colormat = [[0 0 0];[0 0 1];[1 0 0];[nan nan nan]];
if any(IDXNew==0) || any(isnan(IDXNew)) %not sure why this was here.... but here we are
    IDXNew(IDXNew==0 | isnan(IDXNew)) = 4;
end
coloridx = colormat(IDXNew',:);
subplot(1,3,[2,3])
hold all
scatter3(broadbandSlowWave,thratio,EMG,10,coloridx,'filled')
view(133.7,18.8);
grid on
xlabel('Broadband SW');ylabel('Narrowband Theta');zlabel('EMG')
subplot(3,3,1)
hold on
bar(swhistbins,swhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
plot([swthresh swthresh],[0 max(swhist)],'r')
xlabel('Slow Wave')
title('Step 1: Broadband for NREM')
subplot(3,3,4)
hold on
bar(EMGhistbins,EMGhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
plot([EMGthresh EMGthresh],[0 max(EMGhist)],'r')
xlabel('EMG')
title('Step 2: EMG for Muscle Tone')
subplot(3,3,7)
hold on
bar(THhistbins,THhist,'FaceColor','none','barwidth',0.9,'linewidth',2)
plot([THthresh THthresh],[0 max(THhist)],'r')
xlabel('Theta')
title('Step 3: Theta for REM')
        
% saveas(gcf,[figloc,recordingname,'_SSCluster3D'],'jpeg')