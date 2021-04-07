function SpikeParameters = CreateLib_WFParam(mice, samplingrate)
%
% Creates an array of spikes' waveforms and their paramteres to classify
% them later in pyramidal cells and interneurons
%
% INPUT
%
%     mice                number of mice from UMazePAG ERC experiment to
%                         take the neurons from
%     samplingrate        sampling rate of recordings (in Hz)
% 
%  OUTPUT
%
%     SpikeParameters     matrix of N * 37 where N is number of waveforms.
%                         Last 32 of 37 columns correspond to 32 points of a waveform
%                         Column 1 is reserved for classification, 2nd
%                         column is FR, 3rd column is half amplitude
%                         duration, 4th column is assymetry index and 5th
%                         column is trough to peak time
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 05/04/2021
% github.com/bryzgalovdm

%% Parameters
resample_factor = 300;

%% Allocate
TotalEpoch = cell(length(mice), 1);
id_ripples = cell(length(mice), 1);

%% Load data
% Get paths
Dir = PathForExperimentsERC('UMazePAG');
Dir = RestrictPathForExperiment(Dir,'nMice',mice);

for imouse = 1:length(Dir.path)
    wave{imouse} = load([Dir.path{imouse}{1} 'MeanWaveform.mat']);
    spikes{imouse} = load([Dir.path{imouse}{1} 'SpikeData.mat'], 'S', 'RippleGroups', 'TT');
end

%% Calculate epoch
for imouse = 1:length(Dir.path)
    lfp = load([Dir.path{imouse}{1} '/LFPData/LFP0.mat']);
    time = Range(lfp.LFP);
    TotalEpoch{imouse} = intervalSet(time(1), time(end));
end
clear lfp

%% Find ripples channel
for imouse = 1:length(Dir.path)
    id_ripples{imouse} = PyramLayerSorting(spikes{imouse}.S, spikes{imouse}.RippleGroups, spikes{imouse}.TT);
end

%% Calculate firing rate
FR = nan(2e4, 1);
cnt=1;
for imouse = 1:length(Dir.path)
    for ineuron = 1:length(spikes{imouse}.S)
        if sum(ineuron==id_ripples{imouse}) > 0
            FR(cnt) = GetFiringRate(spikes{imouse}.S(ineuron), TotalEpoch{imouse});
            cnt = cnt+1;
        end
    end
end
FR(isnan(FR)) = [];

%% All waveforms together
cnt=1;
AllWF = nan(32, 2e4);
for imouse = 1:length(Dir.path)
    for ineuron = 1:length(wave{imouse}.W)
        if sum(ineuron==id_ripples{imouse}) > 0
            [~, ChanToUse] = min(min(wave{imouse}.W{ineuron}, [], 2));
            AllWF(:,cnt) = wave{imouse}.W{ineuron}(ChanToUse, :);
            cnt=cnt+1;
        end
    end
end
id_nan = isnan(AllWF(14,:));
AllWF(:,id_nan) = [];

%%
for iwave = 1:size(AllWF, 2)
    WaveToUse = AllWF(:, iwave);
    
    if not(sum(isnan(WaveToUse)) == length(WaveToUse))
        % resample to higher frequency
        WaveToUseResample = resample(WaveToUse, resample_factor, 1);
        
        % normalize amplitude
        WaveToUseResample = WaveToUseResample./(max(WaveToUseResample)-min(WaveToUseResample));
        
        % Trough To Peak
        [valMin,indMin] = min(WaveToUseResample); % find trough
        [~,indPeak] = max(WaveToUseResample(indMin:end)); % find next peak
        WFInfo.TroughToPeakTime(iwave) = (indPeak/samplingrate)/resample_factor;
        
        % Half amplitude duration
        HalfAmp = valMin/2;
        TimeAtHlafAmp(1) = find(WaveToUseResample<HalfAmp,1,'first');
        TimeAtHlafAmp(2) = find(WaveToUseResample<HalfAmp,1,'last');
        WFInfo.HalfAmpDur(iwave) = (TimeAtHlafAmp(2)-TimeAtHlafAmp(1))*5e-5/resample_factor;
        
        % Half width
        DD = diff(WaveToUseResample);
        diffpeak = find(DD(indMin:end) == max(DD(indMin:end)))+indMin;
        DD = DD(diffpeak:end);
        IndMax = find(DD<max(abs(diff(WaveToUseResample)))*0.01,1,'first')+diffpeak;
        if isempty(IndMax)
            IndMax = find(DD<max(abs(diff(WaveToUseResample)))*0.05,1,'first')+diffpeak;
        end
        if WaveToUseResample(IndMax)<0
            if not(isempty(find(WaveToUseResample(IndMax:end)>0,1,'first')+IndMax))
                IndMax = find(WaveToUseResample(IndMax:end)>0,1,'first')+IndMax ;
            end
        end
        
        if not(isempty(IndMax)) && not(isempty(IndMax))
            WFInfo.HalfWidth(iwave) = ((IndMax-indMin)/samplingrate)/resample_factor;
        else
            WFInfo.HalfWidth(iwave) = NaN;
        end
        
        % Area under curve
        WaveToUseResampleTemp = WaveToUseResample(indMin:end);
        valzero = find(WaveToUseResampleTemp>0,1,'first');
        WaveToCalc = WaveToUseResampleTemp(valzero:end);
        WFInfo.AreaUnderCurve(iwave) = sum(abs(WaveToCalc));
        if ~isempty(valzero)
            WFInfo.AreaUnderCurveNorm(iwave) = sum(abs(WaveToCalc))./(length(WaveToUseResample)-valzero);
        else
            WFInfo.AreaUnderCurveNorm(iwave) = 0;
        end
        
        % Assymetry
        MaxBef = max(WaveToUseResample(1:indMin));
        MaxAft = max(WaveToUseResample(indMin:end));
        WFInfo.Assymetry(iwave) = (MaxAft-MaxBef)./(MaxAft+MaxBef);
        
        % The wave itself
        WFInfo.Wave(iwave,:) = WaveToUse;
        
    else
        WFInfo.TroughToPeakTime(iwave) = NaN;
        WFInfo.HalfAmpDur(iwave) = NaN;
        WFInfo.HalfWidth(iwave) = NaN;
        WFInfo.AreaUnderCurve(iwave) = NaN;
        WFInfo.AreaUnderCurveNorm(iwave)=NaN;
        WFInfo.Assymetry(iwave)=NaN;
        WFInfo.OldOrNew(iwave)=(iwave<=size(LibraryWF,2));
    end
end
    
%% Output
SpikeParameters = [FR (WFInfo.HalfAmpDur/range(WFInfo.HalfAmpDur))' ...
    WFInfo.Assymetry' (WFInfo.TroughToPeakTime/range(WFInfo.TroughToPeakTime))'];


%% Save Output
AllParams = [nan(size(SpikeParameters, 1), 1) SpikeParameters WFInfo.Wave];

save([dropbox filesep 'Kteam' filesep 'PrgMatlab' filesep 'WaveFormLibraryHip.mat'], 'AllParams');
end