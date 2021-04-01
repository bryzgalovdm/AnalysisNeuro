%CalcBasicNeuronInfo - Calculate basic properties of sorted clusters.
%
% Gets indexes of MUA and SUA, computes firing rate, classify neurons into
% putative Int/Pyr, computes isolation distances , and fetches subjective
% quality from .clu files
%
%  USAGE
%
%    [MatInfoNeurons, BasicNeuronInfo] = CalcBasicNeuronInfo(Dir, ploto)
%
%    Dir       Directoty to find files
%    ploto     if 1 will plot the properties and save the figure; 0 - not
%
%  OUTPUT
%
%    MatInfoNeurons   A matrix with number of neurons (:,1), firing rate
%    (:,2) and putative class (:,3) (Pyr == 1, Int == -1, Amb = +-0.5)
%    BasicNeuronInfo    Structure with all the Basic Info
%
%       See
%
%       MakeNeuronInfoData_KJ, PlotBasicSpikeData, GetFiringRate,
%       MakeData_ClassifySpikeWaveforms, IsolationDistance
%
% Copyright (C) 2018 by Dmitri Bryzgalov
% Changed by Dmitri Bryzgalov 28.03.2019
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

%% TODO:  - divide power for OB (FindStrongOsc)


function [MatInfoNeurons, BasicNeuronInfo] = CalcBasicNeuronInfo(Dir, ploto)
%% Parameters

% Folders with data
if ~exist('Dir', 'var')
    Dir = '/media/nas5/ProjetERC2/Mouse-994/20191013/PagExp/_Concatenated/';
end
% Dir = PathForExperimentsERC_Dima('AllSpikes');

if ~exist('ploto', 'var')
    ploto = 1;
end

%%

disp(' ')
disp('****************************************************************')
cd(Dir)
disp(pwd)

%% load data
load('SpikeData.mat');
load('MeanWaveform');
load('ExpeInfo.mat');

% Response to ripples
if exist('SWR.mat', 'file') == 2
    load('SWR.mat');
else
    load('Ripples.mat');
end

% Different states
if exist('SleepScoring_OBGamma.mat', 'file') == 2
    load('SleepScoring_OBGamma.mat','Wake', 'SWSEpoch', 'REMEpoch');
elseif exist('SleepScoring_Accelero.mat', 'file') == 2
    load('SleepScoring_Accelero.mat','Wake', 'SWSEpoch', 'REMEpoch');
else
    error('Sleep score your data first')
end

% Modulation theta
if exist('ChannelsToAnalyse/dHPC_deep.mat','file')==2
    load('ChannelsToAnalyse/dHPC_deep.mat')
    channel_hpc=channel;
elseif exist('ChannelsToAnalyse/dHPC_rip.mat','file')==2
    load('ChannelsToAnalyse/dHPC_rip.mat')
    channel_hpc=channel;
elseif exist('ChannelsToAnalyse/dHPC_sup.mat','file')==2
    load('ChannelsToAnalyse/dHPC_sup.mat')
    channel_hpc=channel;
else
    disp('No HPC channel, cannot do theta modulation');
end

LFP_hpc = load(['LFPData/LFP' num2str(channel_hpc) '.mat']);
FilLFP_hpc =FilterLFP(LFP_hpc.LFP,[5 10],1024);
time = Range(FilLFP_hpc);

zrTheta = hilbert(Data(FilLFP_hpc));
power=abs(zrTheta);
powerTsd=tsd(Range(FilLFP_hpc),power);
th=percentile(power,50);

ThetaEpoch = thresholdIntervals(Restrict(powerTsd,Wake), th);

% Modulation OB
if exist('ChannelsToAnalyse/Bulb_deep.mat','file')==2
    load('ChannelsToAnalyse/Bulb_deep.mat')
    channel_ob=channel;
end

if exist('channel_ob','var')
    LFP_ob = load(['LFPData/LFP' num2str(channel_ob) '.mat']);
    FilLFP_ob =FilterLFP(LFP_ob.LFP,[3 6],1024);
    
    
    zr4Hz = hilbert(Data(FilLFP_ob));
    power=abs(zr4Hz);
    powerTsd=tsd(Range(FilLFP_ob),power);
    FourHz=percentile(power,75);
    
    FourHzEpoch = thresholdIntervals(Restrict(powerTsd,Wake), FourHz);
end

% Filename
flnme = ['M' num2str(ExpeInfo.nmouse) '_' num2str(ExpeInfo.date) '_' ExpeInfo.SessionType];

%% Get indexes of MUA and SUA
for i = 1:length(TT)
    if TT{i}(2) == 1
        idx_MUA(i) = i;
    else
        idx_SUA(i) = i;
    end
end
idx_MUA(idx_MUA==0) = [];
idx_SUA(idx_SUA==0) = [];

%% Get info

% Calculate overall firing rate and during different states
for j=1:length(S)
    firingrates(j) = length(Range(S{j}))/(time(end)-time(1))*1e4;
    
    FRWake(j) = length(Range(Restrict(S{j},Wake)))/sum((End(Wake)-Start(Wake)))*1e4;
    FR_SWS(j) = length(Range(Restrict(S{j},SWSEpoch)))/sum((End(SWSEpoch)-Start(SWSEpoch)))*1e4;
    FR_REM(j) = length(Range(Restrict(S{j},REMEpoch)))/sum((End(REMEpoch)-Start(REMEpoch)))*1e4;
    FR_WakeTheta(j) = length(Range(Restrict(S{j},and(Wake, ThetaEpoch))))/sum((End(and(Wake, ThetaEpoch))-...
        Start(and(Wake, ThetaEpoch))))*1e4;
    FR_WakeNoTheta(j) = length(Range(Restrict(S{j},(Wake-and(Wake, ThetaEpoch)))))/sum((End((Wake-and(Wake, ThetaEpoch)))-...
        Start((Wake-and(Wake, ThetaEpoch)))))*1e4;
end
% Classify neurons in INT and PYR
if exist('NeuronClassification.mat')~=2
    [UnitID,AllParamsNew,WFInfo,BestElec,figid] = MakeData_ClassifySpikeWaveforms(W,'/home/mobsrick/Dropbox/Kteam/',1);
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    saveas(figid,'NeuronClasses.fig');
    saveFigure(figid,'NeuronClasses',pwd);
    close(figid);
else
    load('NeuronClassification.mat', 'UnitID');
end

% Calculate theta modulation for each unit
for i=1:length(S)
    %[ph{i},mu(i), Kappa(i),pval(i)]=ModulationTheta(S{i},FilLFP_hpc,ThetaEpoch,10,1); % old version of the code
    [ph_temp{i},mu_temp{i},Kappa_temp{i},pval_temp{i},Z_temp{i}]=UnitModulationLFP_SB(S{i},FilLFP_hpc,ThetaEpoch,10,0,1);
end
% Make the output easy to use
for i=1:length(S)
    ph.Nontransf{i} = ph_temp{i}.Nontransf;
    ph.Transf{i} = ph_temp{i}.Transf;
    
    mu.Nontransf(i) = mu_temp{i}.Nontransf;
    mu.Transf(i) = mu_temp{i}.Transf;
    
    Kappa.Nontransf(i) = Kappa_temp{i}.Nontransf;
    Kappa.Transf(i) = Kappa_temp{i}.Transf;
    
    pval.Nontransf(i) = pval_temp{i}.Nontransf;
    pval.Transf(i) = pval_temp{i}.Transf;
    
    Z.Nontransf(i) = Z_temp{i}.Nontransf;
    Z.Transf(i) = Z_temp{i}.Transf;
end
clear ph_temp mu_temp Kappa_temp pval_temp Z_temp

if exist('channel_ob','var')
    % Calculate 4Hz modulation for each unit
    for i=1:length(S)
        %[ph{i},mu(i), Kappa(i),pval(i)]=ModulationTheta(S{i},FilLFP_hpc,ThetaEpoch,10,1); % old version of the code
        [ph_temp{i},mu_temp{i},Kappa_temp{i},pval_temp{i},Z_temp{i}]=UnitModulationLFP_SB(S{i},FilLFP_ob,FourHzEpoch,10,0,1);
    end
    % Make the output easy to use
    for i=1:length(S)
        ph4Hz.Nontransf{i} = ph_temp{i}.Nontransf;
        ph4Hz.Transf{i} = ph_temp{i}.Transf;
        
        if ~isempty(mu_temp{i}.Nontransf) || ~isempty(mu_temp{i}.Transf)
            mu4Hz.Nontransf(i) = mu_temp{i}.Nontransf;
            mu4Hz.Transf(i) = mu_temp{i}.Transf;
        else
            mu4Hz.Nontransf(i) = NaN;
            mu4Hz.Transf(i) = NaN;
        end
        
        if ~isempty(Kappa_temp{i}.Nontransf) || ~isempty(Kappa_temp{i}.Transf)
            Kappa4Hz.Nontransf(i) = Kappa_temp{i}.Nontransf;
            Kappa4Hz.Transf(i) = Kappa_temp{i}.Transf;
        else
            Kappa4Hz.Nontransf(i) = NaN;
            Kappa4Hz.Transf(i) = NaN;
        end
        
        if ~isempty(pval_temp{i}.Nontransf) || ~isempty(pval_temp{i}.Transf)
            pval4Hz.Nontransf(i) = pval_temp{i}.Nontransf;
            pval4Hz.Transf(i) = pval_temp{i}.Transf;
        else
            pval4Hz.Nontransf(i) = NaN;
            pval4Hz.Transf(i) = NaN;
        end
        
        if ~isempty(Z_temp{i}.Nontransf) || ~isempty(Z_temp{i}.Transf)
            Z4Hz.Nontransf(i) = Z_temp{i}.Nontransf;
            Z4Hz.Transf(i) = Z_temp{i}.Transf;
        else
            Z4Hz.Nontransf(i) = NaN;
            Z4Hz.Transf(i) = NaN;
        end
    end
    
    clear ph_temp mu_temp Kappa_temp pval_temp Z_temp
end

% Calculate response to ripples
for i=1:length(S)
    [CC_neurip(i,:),BT_neurip(i,:)]=CrossCorr(ripples(:,2)*1e4,Range(S{i}),1,300);
end


%% save

MatInfoNeurons(:,1) = 1:length(S);
MatInfoNeurons(:,2) = firingrates';
MatInfoNeurons(:,3) = UnitID(:,1); % Pyr == 1, Int == -1, Amb = +-0.5

BasicNeuronInfo.idx_MUA = idx_MUA;
BasicNeuronInfo.idx_SUA = idx_SUA;
BasicNeuronInfo.firingrate = firingrates;
BasicNeuronInfo.FRWake = FRWake;
BasicNeuronInfo.FR_SWS = FR_SWS;
BasicNeuronInfo.FR_REM = FR_REM;
BasicNeuronInfo.FR_WakeTheta = FR_WakeTheta;
BasicNeuronInfo.FR_WakeNoTheta = FR_WakeNoTheta;
BasicNeuronInfo.neuroclass = UnitID(:,1);
BasicNeuronInfo.phasestheta = ph;
BasicNeuronInfo.mutheta = mu;
BasicNeuronInfo.kappatheta = Kappa;
BasicNeuronInfo.pvaltheta = pval;
BasicNeuronInfo.Z = Z;
if exist('channel_ob','var')
    BasicNeuronInfo.phases4Hz = ph4Hz;
    BasicNeuronInfo.mu4Hz = mu4Hz;
    BasicNeuronInfo.kappa4Hz = Kappa4Hz;
    BasicNeuronInfo.pval4Hz = pval4Hz;
    BasicNeuronInfo.Z4Hz = Z4Hz;
end
BasicNeuronInfo.CC_neurip = CC_neurip;
BasicNeuronInfo.BT_neurip = BT_neurip;
BasicNeuronInfo.ThetaEpoch = ThetaEpoch;

if exist('SpikeData.mat')==2
    save('SpikeData', 'MatInfoNeurons', 'BasicNeuronInfo','-append');
else
    save('SpikeData', 'MatInfoNeurons', 'BasicNeuronInfo');
end

%% Plot a figure
if ploto
    PlotBasicSpikeData(BasicNeuronInfo,Quality,1)
end

end