function [UnitID,AllParamsNew,WFInfo,BestElec,figid] = ClassifyHippocampalWaveforms(W, S, LFP, DropBoxLocation,PlotOrNot, varargin)
%
% Creates an array of spikes' waveforms and their paramteres to classify
% them later in pyramidal cells and interneurons
%
% INPUT
%
%     W                   cell array with mean waveforms to be classified
%     S                   tsdArray with all spike time of these waveforms
%     LFP                 any LFP trace corresponding to this recording
%     DropBoxLocation     path to the folder where library waveforms are
%                         stored
%     PlotOrNot           boolean whether to plot classification figure
%     recompute           boolean whether to recompute existing
%                         classification
%     addlib              boolean whether to add newly classified
%                         waverforms to the library
% 
%  OUTPUT
%
%     UnitID              N * 2 matrix where N is number of newly classified neurons,
%                         1st column is classification id (>0 - pyr, <0 - int), 2nd
%                         column is distance to average class member
%     AllParamsNew        N * 4 matrix where N is number of newly classified neurons,
%                         1st column is FR, 2nd column is half amplitude
%                         duration, 3rd column is assymetry index and 4th
%                         column is trough to peak time
%     WFInfo              structure with more parameters of waveforms
%     BestElec            number of best electrode for every waveforms
%     figid               figure id (empty if PlotOrNot=false)
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 05/04/2021
% github.com/bryzgalovdm

%% INITIATION

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property.']);
    end
    switch(lower(varargin{i}))
        case 'recompute'
            recompute = varargin{i+1};
            if recompute~=0 && recompute ~=1
                error('Incorrect value for property ''recompute''.');
            end
        case 'addlib'
            addlib = varargin{i+1};
            if addlib~=0 && addlib ~=1
                error('Incorrect value for property ''addlib''.');
            end
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) '''.']);
    end
end

%check if exist and assign default value if not
if ~exist('recompute','var')
    recompute=0;
end
if ~exist('PlotOrNot','var')
    PlotOrNot = 0;
end
if ~exist('addlib','var')
    addlib = 0;
end

%check if already exist
if ~recompute
    if exist('NeuronClassification.mat','file')==2
        disp('Already computed! ')
        return
    end
end

% Parameters
resample_factor = 300;

% Get waveforms from electrods with max amplitude
for ww=1:length(W)
    clear Peak
    for elec=1:4
        try
            Peak{ww}(elec) = min(W{ww}(elec,:));
        end
    end
    [~,BestElec{ww}] = min(Peak{ww});
    NewWF(ww,:) = W{ww}(BestElec{ww},:);
end
NewWF=NewWF';

% Load the  library of waveforms
load(fullfile(DropBoxLocation,'PrgMatlab', 'WaveFormLibraryHip.mat'))
LibraryWF = (AllParams(:,end-31:end)');

% Concatenate new and all WF
AllWF = [LibraryWF,NewWF];

%% Calculate firing rate
time = Range(LFP);
TotalEpoch = intervalSet(time(1), time(end));
FR = GetFiringRate(S, TotalEpoch);
AllFR = [AllParams(:,2); FR];

%% DEFINE WAVEFORM PARAMETERS
for k = 1:size(AllWF,2)
    
    if not(sum(isnan(AllWF(:,k))) == length(AllWF(:,k)))
        % resample to higher frequency
        WaveToUseResample(k,:) = resample(AllWF(:,k),resample_factor,1);
        
        % normalize amplitude
        WaveToUseResample(k,:) = WaveToUseResample(k,:)./(max(WaveToUseResample(k,:))-min(WaveToUseResample(k,:)));
        
        % Trough To Peak
        [valMin,indMin] = min(WaveToUseResample(k,:)); % find trough
        [~,indPeak] = max(WaveToUseResample(k,indMin:end)); % find next peak
        WFInfo.TroughToPeakTime(k) = indPeak*5e-5/resample_factor;
        
        % Half amplitude duration
        HalfAmp = valMin/2;
        TimeAtHlafAmp(1) = find(WaveToUseResample(k,:)<HalfAmp,1,'first');
        TimeAtHlafAmp(2) = find(WaveToUseResample(k,:)<HalfAmp,1,'last');
        WFInfo.HalfAmpDur(k) = (TimeAtHlafAmp(2)-TimeAtHlafAmp(1))*5e-5/resample_factor;
        
        % Half width
        DD = diff(WaveToUseResample(k,:));
        diffpeak = find(DD(indMin:end) == max(DD(indMin:end)))+indMin;
        DD = DD(diffpeak:end);
        IndMax = find(DD<max(abs(diff(WaveToUseResample(k,:))))*0.01,1,'first')+diffpeak;
        if isempty(IndMax)
            IndMax = find(DD<max(abs(diff(WaveToUseResample(k,:))))*0.05,1,'first')+diffpeak;
        end
        if WaveToUseResample(k,IndMax)<0
            if not(isempty(find(WaveToUseResample(k,IndMax:end)>0,1,'first')+IndMax))
                IndMax = find(WaveToUseResample(k,IndMax:end)>0,1,'first')+IndMax ;
            end
        end
        
        if not(isempty(IndMax)) && not(isempty(IndMax))
            WFInfo.HalfWidth(k) = ((IndMax-indMin)*5e-5)/resample_factor;
        else
            WFInfo.HalfWidth(k) = NaN;
        end
        
        % Area under curve
        WaveToUseResampleTemp = WaveToUseResample(k,indMin:end);
        valzero = find(WaveToUseResampleTemp>0,1,'first');
        WaveToCalc = WaveToUseResampleTemp(valzero:end);
        WFInfo.AreaUnderCurve(k) = sum(abs(WaveToCalc));
        if ~isempty(valzero)
            WFInfo.AreaUnderCurveNorm(k) = sum(abs(WaveToCalc))./(length(WaveToUseResample(k,:))-valzero);
        else
            WFInfo.AreaUnderCurveNorm(k) = 0;
        end
        
        % Assymetry
        MaxBef = max(WaveToUseResample(k,1:indMin));
        MaxAft = max(WaveToUseResample(k,indMin:end));
        WFInfo.Assymetry(k) = (MaxAft-MaxBef)./(MaxAft+MaxBef);
        
        % is this for a new (0) or a library neuron (1)
        WFInfo.OldOrNew(k)=(k<=size(LibraryWF,2));
    else
        WFInfo.TroughToPeakTime(k) = NaN;
        WFInfo.HalfAmpDur(k) = NaN;
        WFInfo.HalfWidth(k) = NaN;
        WFInfo.AreaUnderCurve(k) = NaN;
        WFInfo.AreaUnderCurveNorm(k)=NaN;
        WFInfo.Assymetry(k)=NaN;
        WFInfo.OldOrNew(k)=(k<=size(LibraryWF,2));
    end
end


%% CLASSIFICATION OF UNITS
% Cluster
rmpath(fullfile(DropBoxLocation,'PrgMatlab/Fra/UtilsStats'));
Params = [AllFR (WFInfo.HalfAmpDur/range(WFInfo.HalfAmpDur))' ...
    WFInfo.Assymetry' (WFInfo.TroughToPeakTime/range(WFInfo.TroughToPeakTime))'];
NeuronClassif = kmeans(Params,2);
addpath(fullfile(DropBoxLocation,'PrgMatlab/Fra/UtilsStats'));

% 1 is Pyr, -1 is Int
NumofOnes = sum(NeuronClassif == 1)/length(NeuronClassif);
if NumofOnes>0.5
    NeuronClassif(NeuronClassif == 1) = 1;
    NeuronClassif(NeuronClassif == 2) = -1;
else
    NeuronClassif(NeuronClassif == 1) = -1;
    NeuronClassif(NeuronClassif == 2) = 1;
end


% Look at distance to average WF to find the WF that are hard to classify
% Pyramidal neurons
MeanPyr = mean(AllWF(:,NeuronClassif == 1)')';
PyrNeurons = find(NeuronClassif == 1);
for ff = 1:length(PyrNeurons)
    DistToMeanPyr(ff) = sum(abs(AllWF(:,PyrNeurons(ff))-MeanPyr));
end
LimPyr = mean(DistToMeanPyr)+2*std(DistToMeanPyr);
AmbigPyr = PyrNeurons((DistToMeanPyr>LimPyr));
NeuronClassif(AmbigPyr) = 0.5;
NeuronClassif(PyrNeurons,2)=DistToMeanPyr;

MeanInt = mean(AllWF(:,NeuronClassif == -1)')';
IntNeurons = find(NeuronClassif == -1);
for ff = 1:length(IntNeurons)
    DistToMeanInt(ff) = sum(abs(AllWF(:,IntNeurons(ff))-MeanInt));
end
LimInt = mean(DistToMeanInt)+2*std(DistToMeanInt);
AmbigInt = IntNeurons((DistToMeanInt>LimInt));
NeuronClassif(AmbigInt) = -0.5;
NeuronClassif(IntNeurons,2) = DistToMeanInt;

%Separate the old from the new

UnitIDNew = NeuronClassif(end-size(NewWF,2)+1:end,:);
AllParamsNew = Params(end-size(NewWF,2)+1:end,:);
DatNew = AllWF(:,end-size(NewWF,2)+1:end);

UnitIDOld = NeuronClassif(1:size(LibraryWF,2),:);
AllParamsOld = Params(1:size(LibraryWF,2),:);
DatOld = AllWF(:,1:size(LibraryWF,2));

%result
UnitID = UnitIDNew;


%% SAVE
save NeuronClassification UnitID AllParamsNew BestElec
if addlib
    POld = [UnitIDOld(:,1) AllParamsOld DatOld'];
    PNew = [UnitIDNew(:,1) AllParamsNew DatNew'];
    AllParams = [POld; PNew];
    save([dropbox filesep 'Kteam' filesep 'PrgMatlab' filesep 'WaveFormLibraryHip.mat'], 'AllParams');
end

%% PLOT FIGURE
if PlotOrNot
    figid=figure;
    
    subplot(2,3,[1,2])
    %plot old units
    plot3(AllParamsOld(UnitIDOld == -1,3),AllParamsOld(UnitIDOld == -1,2),AllParamsOld(UnitIDOld == -1,1),'r.','MarkerSize',5)
    hold on,plot3(AllParamsOld(UnitIDOld == 1,3),AllParamsOld(UnitIDOld == 1,2),AllParamsOld(UnitIDOld == 1,1),'b.','MarkerSize',5)
    plot3(AllParamsOld(UnitIDOld == -0.5,3),AllParamsOld(UnitIDOld == -0.5,2),AllParamsOld(UnitIDOld == -0.5,1),'m.','MarkerSize',5)
    plot3(AllParamsOld(UnitIDOld == 0.5,3),AllParamsOld(UnitIDOld == 0.5,2),AllParamsOld(UnitIDOld == 0.5,1),'c.','MarkerSize',5)
    
    % plot new units
    plot3(AllParamsNew(UnitIDNew == -1,3),AllParamsNew(UnitIDNew == -1,2),AllParamsNew(UnitIDNew == -1,1),'ro','MarkerSize',5)
    hold on,plot3(AllParamsNew(UnitIDNew == 1,3),AllParamsNew(UnitIDNew == 1,2),AllParamsNew(UnitIDNew == 1,1),'bo','MarkerSize',5)
    plot3(AllParamsNew(UnitIDNew == -0.5,3),AllParamsNew(UnitIDNew == -0.5,2),AllParamsNew(UnitIDNew == -0.5,1),'mo','MarkerSize',5)
    plot3(AllParamsNew(UnitIDNew == 0.5,3),AllParamsNew(UnitIDNew == 0.5,2),AllParamsNew(UnitIDNew == 0.5,1),'co','MarkerSize',5)
    legend({'Int','PN','IntAmbig','PyrAmbig','Int new','PN new','IntAmbig new','PyrAmbig new'})
    zlabel('Firing rate (Hz)');
    ylabel('Time between halfAmp');
    xlabel('Asymetry');
    
    subplot(2,3,4)
    hist(UnitIDOld(UnitIDOld(:,1)>0,2),50,'k'),
    h  =  findobj(gca,'Type','patch');
    set(h,'FaceColor','b','EdgeColor','w')
    hold on, line([LimPyr LimPyr],get(gca,'ylim'),'color','k','linewidth',3)
    plot(UnitIDNew(UnitIDNew(:,1)>0,2),UnitIDNew(UnitIDNew(:,1)>0,2)*0+max(get(gca,'ylim'))/2,'k*')
    xlabel('black bar  =  mean + 2SD')
    ylabel('Dist To Mean Pyr')
    
    
    subplot(2,3,5)
    hist(UnitIDOld(UnitIDOld(:,1)<0,2),50,'k'),
    h  =  findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    hold on, line([LimInt LimInt],get(gca,'ylim'),'color','k','linewidth',3)
    plot(UnitIDNew(UnitIDNew(:,1)<0,2),UnitIDNew(UnitIDNew(:,1)<0,2)*0+max(get(gca,'ylim'))/2,'k*')
    ylabel('Dist To Mean Int')
    
    subplot(2,3,[3])
    plot(DatOld(:,UnitIDOld(:,1) == 1),'b'),
    hold on
    if ~isempty(PyrNeurons)
        plot(DatNew(:,UnitIDNew(:,1) == 1),'c'),
        plot(DatNew(:,UnitIDNew(:,1) == 0.5),'k'),
    end
    xlabel('blue: old   cyan: new   black: ambig')
    
    subplot(2,3,[6])
    plot(DatOld(:,UnitIDOld(:,1) == -1),'r'),
    hold on
    if ~isempty(IntNeurons)
        plot(DatNew(:,UnitIDNew == -1),'m'),
        plot(DatNew(:,UnitIDNew == -0.5),'k'),
    end
    xlabel('red: old   magenta: new   black: ambig')
else
    figid=[];
end

disp(['Proportion of interneurons new data : ' num2str(length(find(UnitIDNew(:,1) == -1))./(length(find(UnitIDNew(:,1) == -1))+length(find(UnitIDNew(:,1) == 1)))*100) '%'])
disp(['Proportion of interneurons library data: ' num2str(length(find(UnitIDOld(:,1) == -1))./(length(find(UnitIDOld(:,1) == -1))+length(find(UnitIDOld(:,1) == 1)))*100) '%'])


end