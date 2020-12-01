function [bootstrapped, fh] = PreTestBootstrap(NewMouseData, varargin)
%
% This function bootstraps occupancies of shock zone and safe zone for Pre-Tests:
%   It plots the bootstrapping distributions and it's 95% confidence intervals 
%   and places selected mouse in the distribution
%
% INPUT
%
%     NewMouseData   2-elements vector with mean occupancy of ShockZone
%                    (NewMouseData(1)) and SafeZone(NewMouseData(2)) in
%                    PreTests (in ratios!)
%     Nboot          (optional) number of bootstrapps to do (default=1000)
%     ShowFigure     (optional) 1 if you want to plot the figure (default=1)
% 
%  OUTPUT
%
%     boostrapped   matrix (Nboot*2) of bootrapped means for Shock and Safe
%                   zones
%     fh            handle to the figure
% 
%           EXAMPLE:
%               [bootstrapped, fh] = PreTestBootstrap([0.18 0.25]);
%               bootstrapped = PreTestBootstrap([0.18 0.25], 'NBoot', '2000');
%               [bootstrapped, fh] = PreTestBootstrap([0.18 0.25], 'NBoot', '2000', 'ExpType', 'MFB');
%       
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 01/12/2020
% github.com/bryzgalovdm

%% Hyperparameters
nMicePAG = [797 798 828 861 882 905 906 911 912 977 994];
nMiceMFB = [934 941 863 913];

%% Default values of optional arguments
Nboot = 1000;
Type = 'PAG';

%% Optional parameters handling
for i=1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'nboot'
            NBoot = varargin{i+1};
            if ~isa(NBoot, 'double') && NBoot/NBoot ~= 1
                error('Incorrect value for property ''NBoot'' (type ''help PreTestBootstrap'' for details).');
            end
        case 'exptype'
            Type = varargin{i+1};
            if ~strcmp(Type, 'PAG') && ~strcmp(Type, 'MFB')
                error('Incorrect value for property ''ExpType'' (type ''help PreTestBootstrap'' for details).');
            end
    end
end

%% Get data
if strcmp(Type, 'PAG')
    Dir = PathForExperimentsERC_Dima('UMazePAG');
    Dir = RestrictPathForExperiment(Dir,'nMice', nMicePAG);
else
    Dir = PathForExperimentsERC_SL('StimMFBWake');
    Dir = RestrictPathForExperiment(Dir,'nMice', nMiceMFB);
end

a = cell(length(Dir.path),1);
for i = 1:length(Dir.path)
    % PreTests
    a{i} = load([Dir.path{i}{1} '/behavResources.mat'], 'behavResources');
end


%% Find indices of PreTests and PostTest session in the structure
id_Pre = cell(1,length(a));

for i=1:length(a)
    id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
end

%% Calculate average occupancy
% Calculate occupancy de novo
occup = nan(length(Dir.path), 4, 2); % 4 tests, Shock and Safe

for i=1:length(a)
    for k=1:length(id_Pre{i})
        temp = CalculateZoneOccupancy(a{i}.behavResources(id_Pre{i}(k)));
        occup(i,k,1) = temp(1);
        occup(i,k,2) = temp(2);
    end
end

occup_mean = nan(length(Dir.path), 2);
occup_std = nan(length(Dir.path), 2);

for izone = 1:2 % 1 codes for ShockZone, 2 for SafeZone
    occup_mean(:,izone) = mean(squeeze(occup(:,:,izone)),2);
    occup_std(:,izone) = std(squeeze(occup(:,:,izone)),0,2);
end

%% Do bootstrapping

bootstrapped = nan(Nboot,2);  % Nboot bootstraps, Shock and Safe
CI_minus = nan(1,2);
CI_plus = nan(1,2);

for izone = 1:2 % 1 codes for ShockZone, 2 for SafeZone
    bootstrapped(:,izone) = bootstrp(Nboot,@mean,occup_mean(:,izone));
    CI_minus(izone) = prctile(bootstrapped(:,izone),2.5);
    CI_plus(izone) = prctile(bootstrapped(:,izone),97.5);
end


%% Plot
if strcmp(Type, 'PAG')
    titles = {'Shock zone occupancy', 'Safe zone occupancy'};
elseif strcmp(Type, 'MFB')
    titles = {'Reward zone occupancy', 'No-reward zone occupancy'};
end
fh = figure('units', 'normalized', 'outerposition', [0.1 0.2 0.65 0.5]);
for izone = 1:2 % 1 codes for ShockZone, 2 for SafeZone
    subplot(1,2,izone)
    histogram(bootstrapped(:,izone)*100);
    l1 = line([CI_minus(izone) CI_minus(izone)]*100, ylim, 'Color','r', 'LineWidth',2);
    line([CI_plus(izone) CI_plus(izone)]*100, ylim, 'Color','r', 'LineWidth',2);
    l3 = line([NewMouseData(izone) NewMouseData(izone)]*100, ylim, 'Color','b', 'LineWidth',2);
    title(titles{izone}, 'FontWeight', 'bold', 'FontSize', 14);
    xlabel('Occupancy percentage', 'FontWeight', 'bold', 'FontSize', 14);
    if izone ==2
        legend([l1,l3], '95% CI', 'Current Mouse', 'FontWeight', 'bold', 'FontSize', 14);
    end
end