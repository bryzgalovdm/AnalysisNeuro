function ReDefineZones_DB(directories)


for idir = 1:length(directories)
    cd(directories{idir});
    
    if exist('cleanbehavResources.mat', 'file') == 2
        data = load('cleanbehavResources.mat');
    else
        data = load('behavResources.mat');
    end
    
    % Figure
    f = figure;
    imagesc(imfuse(double(data.ref), data.Zone{1}, 'blend')), colormap jet, hold on
    plot(Data(data.CleanXtsd)*data.Ratio_IMAonREAL,Data(data.CleanYtsd)*data.Ratio_IMAonREAL,'color',[0.8 0.8 0.8])
    A = regionprops(data.Zone{1},'Centroid');
    plot(A.Centroid(1),A.Centroid(2),'r.','MarkerSize',50)
    plot(A.Centroid(1),A.Centroid(2),'w*','MarkerSize',10)
    title('Define 7 zones')
    
    % Redefine zones
    [Zone,ZoneLabels] = DefineZones(f, data.mask);
    data.Zone = Zone;
    data.ZoneLabels = ZoneLabels;
    [ZoneEpoch, ZoneIndices, Occup] = RecalculateZones(data);
    data.ZoneEpoch = ZoneEpoch;
    data.ZoneIndices = ZoneIndices;
    data.Occup = Occup;
    
    % Save
    if exist('cleanbehavResources.mat', 'file') == 2
        movefile cleanbehavResources.mat cleanbehavResources_old.mat
        save('cleanbehavResources.mat', '-struct', 'data');
    elseif exist('behavResources.mat', 'file') == 2
        movefile behavResources.mat behavResources_old.mat
        save('behavResources.mat', '-struct', 'data');
    end
        
end
end


%% Define Zones function
function [Zone, ZoneLabels] = DefineZones(fighandle, mask)
% Automatically or manually

figtemp=fighandle;

% added farshock and farnoshock zones (SL - 23/07/2019)
title('Shock'); [x1,y1,Shock,x2,y2]=roipoly; Zone{1}=uint8(Shock);plot(x2,y2,'linewidth',2)
title('NoShock');[x1,y1,NoShock,x2,y2]=roipoly(); Zone{2}=uint8(NoShock);plot(x2,y2,'linewidth',2)
title('Centre');[x1,y1,Centre,x2,y2]=roipoly(); Zone{3}=uint8(Centre);plot(x2,y2,'linewidth',2)
title('CentreShock');[x1,y1,CentreShock,x2,y2]=roipoly(); Zone{4}=uint8(CentreShock);plot(x2,y2,'linewidth',2)
title('CentreNoShock');[x1,y1,CentreNoShock,x2,y2]=roipoly(); Zone{5}=uint8(CentreNoShock);plot(x2,y2,'linewidth',2)
title('FarShock');[x1,y1,FarShock,x2,y2]=roipoly(); Zone{8}=uint8(FarShock);plot(x2,y2,'linewidth',2)
title('FarNoShock');[x1,y1,FarNoShock,x2,y2]=roipoly(); Zone{9}=uint8(FarNoShock);plot(x2,y2,'linewidth',2)
stats = regionprops(mask, 'Area');
tempmask=mask;
AimArea=stats.Area/2;
ActArea=stats.Area;
while AimArea<ActArea
    tempmask=imerode(tempmask,strel('disk',1));
    stats = regionprops(tempmask, 'Area');
    ActArea=stats.Area;
end
Zone{6}=uint8(tempmask); % I]
Zone{7}=uint8(mask-tempmask);% Outside area
ZoneLabels={'Shock','NoShock','Centre','CentreShock','CentreNoShock','Inside','Outside','FarShock','FarNoShock'};

close(fighandle)

end


%% Recalculate zone vars
function [ZoneEpoch, ZoneIndices, Occup] = RecalculateZones(data)

Xtemp=Data(data.CleanXtsd); Ytemp=Data(data.CleanYtsd); T1=Range(data.CleanYtsd);
if isfield(data, 'Zone')
    XXX = floor(Data(data.CleanYtsd)*data.Ratio_IMAonREAL);
        XXX(isnan(XXX)) = 240;
        YYY = floor(Data(data.CleanXtsd)*data.Ratio_IMAonREAL);
        YYY(isnan(YYY)) = 320;
    for t=1:length(data.Zone)
        ZoneIndices{t}=find(diag(data.Zone{t}(XXX,YYY)));
        Xtemp2=Xtemp*0;
        Xtemp2(ZoneIndices{t})=1;
        ZoneEpoch{t}=thresholdIntervals(tsd(T1,Xtemp2),0.5,'Direction','Above');
        Occup(t)=size(ZoneIndices{t},1)./size(Data(data.CleanXtsd),1);
        
    end
end

end