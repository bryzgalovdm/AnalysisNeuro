function TrajAroundStim(mice, expe, varargin)
%
% This functions plot trajectories of a selected mouse in pretests,
% condiioning sessions and posttests as well as bar plots of shock-safe
%  zone occupancies
%
%  INPUT
%       
%       Dir             directory with behavResources.mat
%       ISPCDriven      true if the mouse was stimulated using BCI and info
%                           about stims is stored only in StimEpoch (boolean)
%       IsSave          true if you want to save the resulting figure in
%                           Dir (boolean - default = false)
%       NumTest         number of tests to plot the trajectories (default = 4)
% 
% 
%  OUTPUT
%
%       Figure
%
%       See
%   
%       BehaviorERC, PathForExperimentERC_Dima
% 
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 27/10/2020, based on my own script coded in 2019
% github.com/bryzgalovdm


%% Default parameters
IsSave = false;
numtest=4;
IsColorTests = false;

%% Optional parameters
for i=1:2:length(varargin)    
    switch(lower(varargin{i})) 
        case 'issave'
            IsSave = varargin{i+1};
            if IsSave ~= 1 && IsSave ~= 0
                error('Incorrect value for property ''IsSave'' (type ''help ExampleTrajectory_DB'' for details).');
            end
        case 'numtest'
            numtest =  varargin{i+1};
            if ~isa(numtest, 'double')
                error('Incorrect value for property ''NumTest'' (type ''help ExampleTrajectory_DB'' for details).');
            end
        case 'colortests'
            IsColorTests =  varargin{i+1};
            if IsColorTests ~= 1 && IsColorTests ~= 0
                error('Incorrect value for property ''ColorTests'' (type ''help ExampleTrajectory_DB'' for details).');
            end    
    end
end

%% Manage experiment
if strcmp(expe, 'PAG')
    fetchpaths = 'UMazePAG';
    FigName = 'OccupCondPAG';
elseif strcmp(expe, 'MFB')
    fetchpaths = 'StimMFBWake';
    FigName = 'OccupCondMFB';
elseif strcmp(expe, 'Novel')
    fetchpaths = 'Novel';
    FigName = 'OccupNovel';
end

%% Allocate data arrays
% Get paths of each individual mouse
Dir = PathForExperimentsERC(fetchpaths);
Dir = RestrictPathForExperiment(Dir,'nMice',mice);
numsessions = CountNumSesionsERC(Dir);

% Data
a = cell(length(numsessions), 1);

% Auxiliary
Epochs = cell(length(numsessions), 1);
OccupMap = cell(length(numsessions), 1);
OccupMap_large = cell(length(numsessions), 1);

%% load data
cnt=1;
for imouse=1:length(Dir.path)
    for isession = 1:length(Dir.path{imouse})
       a{cnt} = load([Dir.path{imouse}{isession} 'behavResources.mat'],...
            'behavResources', 'StimEpoch', 'SessionEpoch');
        cnt=cnt+1;
    end
end


%% Find necessary tests
id_Pre = cell(1,length(a));
id_Cond = cell(1,length(a));
id_Post = cell(1,length(a));

for i=1:numsessions
    id_Pre{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPre');
    id_Cond{i} = FindSessionID_ERC(a{i}.behavResources, 'Cond');
    id_Post{i} = FindSessionID_ERC(a{i}.behavResources, 'TestPost');
end


% Epochs
for isession = 1:numsessions
    if strcmp(expe, 'PAG')
        for itest = 1:length(id_Cond{isession})
            Epoch{isession}{itest} = intervalSet(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1,1)*1e4,...
                a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1,1)*1e4+4e4);
        end
    elseif strcmp(expe, 'MFB')
        for itest = 1:length(id_Cond{isession})
            Epoch{isession}{itest} = intervalSet(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1,1)*1e4-4e4,...
                a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1,1)*1e4);
        end
    end
end

%% Plot 
fh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for isession = 1:numsessions
    subplot(numsessions, 3, (isession-1)*3+1)
    hold on
    for itest = 1:length(id_Pre{isession})
        plot(Data(a{isession}.behavResources(id_Pre{isession}(itest)).AlignedXtsd),...
            Data(a{isession}.behavResources(id_Pre{isession}(itest)).AlignedYtsd), 'o',...
            'Color', [0.2 0.2 0.2]);
    end
    set(gca,'XTick', [], 'XTickLabel',{},'YTick', [], 'YTickLabel',{});
    ylim([0 1])
    xlim([0 1])
    title('Pre')
    
    subplot(numsessions, 3, (isession-1)*3+2)
    hold on
    for itest = 1:length(id_Cond{isession})
        plot(Data(Restrict(a{isession}.behavResources(id_Cond{isession}(itest)).AlignedXtsd, Epoch{isession}{itest})),...
            Data(Restrict(a{isession}.behavResources(id_Cond{isession}(itest)).AlignedYtsd, Epoch{isession}{itest})), 'o',...
            'Color', [0.2 0.2 0.2]);
        
        if strcmp(expe, 'PAG')
            tempX = Data(a{isession}.behavResources(id_Cond{isession}(itest)).AlignedXtsd);
            tempY = Data(a{isession}.behavResources(id_Cond{isession}(itest)).AlignedYtsd);
            plot(tempX(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1),tempY(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1),...
                'p','Color','r','MarkerFaceColor','red','MarkerSize',4);
            clear tempX tempY
        elseif strcmp(expe, 'MFB')
            tempX = Data(a{isession}.behavResources(id_Cond{isession}(itest)).AlignedXtsd);
            tempY = Data(a{isession}.behavResources(id_Cond{isession}(itest)).AlignedYtsd);
            plot(tempX(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1),tempY(a{isession}.behavResources(id_Cond{isession}(itest)).PosMat(:,4)==1),...
                'p','Color','g','MarkerFaceColor','green','MarkerSize',4);
            clear tempX tempY
        end
    end
    set(gca,'XTick', [], 'XTickLabel',{},'YTick', [], 'YTickLabel',{});
    ylim([0 1])
    xlim([0 1])
    title('Cond')
    
    subplot(numsessions, 3, (isession-1)*3+3)
    hold on
    for itest = 1:length(id_Post{isession})
        plot(Data(a{isession}.behavResources(id_Post{isession}(itest)).AlignedXtsd),...
            Data(a{isession}.behavResources(id_Post{isession}(itest)).AlignedYtsd), 'o'...
            'Color', [0.2 0.2 0.2]);
    end
    set(gca,'XTick', [], 'XTickLabel',{},'YTick', [], 'YTickLabel',{});
    ylim([0 1])
    xlim([0 1])
    title('Post')
    
end

end
