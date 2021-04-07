function ClassifyHipneurons(mice)
%
% Classifies library data and creates the library file that would be used
% as a baseline for future classification
%
% INPUT
%
%     mice                number of mice from UMazePAG ERC experiment to
%                         take the neurons from
% 
%  OUTPUT
%
%     ~                   the file is saved into /PrgMatlab/WaveFormLibraryHip.mat
%
% Coded by Dima Bryzgalov, MOBS team, Paris, France
% 05/04/2021
% github.com/bryzgalovdm


%% Get spike parameters
SpikeParameters = CreateLib_WFParam(mice, 2e4);

%% Do Clustering
rmpath([dropbox filesep '/Kteam/PrgMatlab/Fra/UtilsStats']);
NeuronClassif=kmeans(SpikeParameters(:,1:3),2);
addpath([dropbox filesep '/Kteam/PrgMatlab/Fra/UtilsStats']);
% Identify the PN cell group as the most numerous one
NumofOnes = sum(NeuronClassif == 1)/length(NeuronClassif);
if NumofOnes>0.5
    NeuronClassif(NeuronClassif == 1) = 1;
    NeuronClassif(NeuronClassif == 2) = -1;
else
    NeuronClassif(NeuronClassif == 1) = -1;
    NeuronClassif(NeuronClassif == 2) = 1;
end

%% Load library and put the classification tags to the first column
load([dropbox filesep 'Kteam' filesep 'PrgMatlab' filesep 'WaveFormLibraryHip.mat']);
AllParams(:,1) = NeuronClassif;

%% Look at distance to average WF to find the WF that are hard to classify
AllWF = (AllParams(:,end-31:end))';
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

% Record
AllParams(:,1) = NeuronClassif(:,1);

%% Plot
figid=figure;

subplot(2,3,[1,2])
% Plot old units
plot3(AllParams(NeuronClassif == -1,4),AllParams(NeuronClassif == -1,3),AllParams(NeuronClassif == -1,2),'r.','MarkerSize',5)
hold on,plot3(AllParams(NeuronClassif == 1,4),AllParams(NeuronClassif == 1,3),AllParams(NeuronClassif == 1,2),'b.','MarkerSize',5)
plot3(AllParams(NeuronClassif == -0.5,4),AllParams(NeuronClassif == -0.5,3),AllParams(NeuronClassif == -0.5,2),'m.','MarkerSize',5)
plot3(AllParams(NeuronClassif == 0.5,4),AllParams(NeuronClassif == 0.5,3),AllParams(NeuronClassif == 0.5,2),'c.','MarkerSize',5)
legend({'Int','PN','IntAmbig','PyrAmbig','Int new','PN new','IntAmbig new','PyrAmbig new'})
zlabel('Firing rate (Hz)');
ylabel('Time between halfAmp');
xlabel('Asymetry');

subplot(2,3,4)
hist(NeuronClassif(NeuronClassif(:,1)>0,2),50,'k'),
h  =  findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w')
hold on, line([LimPyr LimPyr],get(gca,'ylim'),'color','k','linewidth',3)
xlabel('black bar  =  mean + 2SD')
ylabel('Dist To Mean Pyr')


subplot(2,3,5)
hist(NeuronClassif(NeuronClassif(:,1)<0,2),50,'k'),
h  =  findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')
hold on, line([LimInt LimInt],get(gca,'ylim'),'color','k','linewidth',3)
ylabel('Dist To Mean Int')

subplot(2,3,[3])
plot(AllWF(:,NeuronClassif(:,1) == 1),'b'),
hold on
if ~isempty(PyrNeurons)
    plot(AllWF(:,NeuronClassif(:,1) == 0.5),'k'),
end
xlabel('blue: old  black: ambig')

subplot(2,3,[6])
plot(AllWF(:,NeuronClassif(:,1) == -1),'r'),
hold on
if ~isempty(IntNeurons)
    plot(AllWF(:,NeuronClassif == -0.5),'k'),
end
xlabel('red: old   magenta: new   black: ambig')

%% Save the data
save([dropbox filesep 'Kteam' filesep 'PrgMatlab' filesep 'WaveFormLibraryHip.mat'], 'AllParams');

end