%
% Summarizes behavioral state and speed data over space and time. Data is
% reduced by binning data in sapce and time according to BinSettings.
%
% USAGE:
%   BinFileList = WormDensity(FileName,BinSettings)
%
%   FileName: single filename of segmented data or cell array of multiple
%               filenames. Select with user input if none given.
%   BinSettings: vector of binsize [time(min) distance(pixels)]
%                default is [0.5 10]

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 14-Apr-2011 09:59:08 
%---------------------------- 

function BinFileList = WormDensity(FileName,BinSettings)

defaultBinSettings = [0.5 10]; % time bin (min), x-y bin (pix)
if nargin < 2 BinSettings = defaultBinSettings; end
if nargin < 1 || isempty(FileName)
    FileName = {};
end
    
%-------- Settings ------ 
SaveData = 1;
s = 4;
segname = '_seg';
BinFileList = {};

%-- Bin Settings ----------------------
if isempty(BinSettings)
    answer = dagetnum({'Time Bin (min):','Area Bin (pixels):'}, ...
                        defaultBinSettings);
    timebin = answer(1).num;  % minutes per bin
    areabin = answer(2).num; % in pixels
end

timebin = BinSettings(1); %0.5  % minutes per bin
areabin = BinSettings(2); %10;  % in pixels
    
clear FileData

if isempty(FileName)
    
    [FileName, PathName] = uigetfile('*_seg.mat','Select segmented track File(s) For Analysis','MultiSelect','on');
    if ~exist('Tracks') && ~iscell(FileName) && FileName(1) == 0 
        errordlg('No file was selected for analysis');
        return;
    end
    cd(PathName);
    if ~iscell(FileName)
        FileName = cellstr(FileName);
    end
end
    
NumFiles = size(FileName,2);
disp(['Batch analyzing ',num2str(NumFiles),' segmented tracks files.']);

    
%-------------------------------------------------------%
%               BIG LOOP FOR EACH FILE                  %
%-------------------------------------------------------%

for FileNum = 1:NumFiles
       
    if exist('PathName')
        FullName = fullfile(PathName, char(FileName(FileNum)));
    else
        FullName = char(FileName(FileNum));
    end

    [PathName,BaseName,Ext] = fileparts(FullName);
    BaseName = BaseName(1:strfind(BaseName,segname)-1); % remove '_seg' from filename
    FileData(FileNum).PathName = PathName;
    FileData(FileNum).BaseName = BaseName;
    FileData(FileNum).RawFile = fullfile(PathName,[BaseName,Ext]);
    FileData(FileNum).SegFile = fullfile(PathName,[BaseName,segname,Ext]);
    
clear('FileInfo','ExpData','TrackData','All*','Animal*','Track','Beh*','*Speed*'); %free memory

%----------------------
% Get experiment data 
%----------------------

disp([datestr(now),': Loading File #',num2str(FileNum),': ',FullName]);

try
    load(FileData(FileNum).SegFile,'Tracks','File','AllData','ExpData');
    disp([datestr(now),'--     ...loaded OK.']);
catch
    disp('Problem Loading File.  Skipping to next File.');
    continue
end

existAll = exist('AllData');

[mp,mf,me] = fileparts(File(ExpData.MovieNum).MovieFile);
FileData(FileNum).MovieFile = fullfile(PathName,[mf,me]);

sizebins = 0:5:max(AllData(:,2));
singleworm = mean(ExpData.TrackSettings(:,3));

if length(singleworm)>0
    singleworm = singleworm(1);
    disp(['Animal size: ',num2str(singleworm),' pixels']);

%     figure(1);clf;bar(sizebins,hist(AllData(:,2),sizebins)); grid on;
%     set(gca,'XTick',0:singleworm:max(sizebins));

    AllData(:,6) = round(AllData(:,2)/singleworm);  % number of worms per centroid
    SkipArena = 0;
else
    disp(['*** NO SINGLE WORMS FOUND ***  Mean object size = ',num2str(mean(AllData(:,2)))]);
    SkipArena = 1;
end

AllData(:,7) = AllData(:,1) / ExpData.FrameRate / 60;  % time in mins

% Proportion of clumping animals
hrange = 0:9;
[nw,wms] = hist(AllData(:,6),hrange);
total = sum(nw.*wms);
nclump = nw.*wms/total;
disp('Percent animals in collisions of n animals');
disp([wms; nclump*100]);

%--- Trim AllData Variable: only 1+ worms/object and correct timing range
prevsize = size(AllData,1);
AllData = AllData(find(AllData(:,6) >= 1),:); % & AllData(:,7) >= StartTime & AllData(:,7) <= EndTime),:);
sizesaving = 1 - size(AllData,1)/prevsize;
disp(['AllData size decreased by: ',num2str(sizesaving*100),'%']);
AllData = single(AllData);

%---set bin parameters
TimeBin = ceil(AllData(:,7)/timebin); tbins = max(TimeBin);
XBin = ceil(AllData(:,3)/areabin); xbins = max(XBin);
YBin = ceil(AllData(:,4)/areabin); ybins = max(YBin);

AllXYT = zeros(xbins, ybins, tbins,'single');
tic;
for i=1:size(AllData,1)
    if XBin(i)>0 & YBin(i)>0
        AllXYT(XBin(i),YBin(i),TimeBin(i)) = AllXYT(XBin(i),YBin(i),TimeBin(i)) + AllData(i,6);
    end
end

All.XYTime = AllXYT;
t = toc; disp(['AllXYT analysis: ',num2str(t),' s']);

All.XY = squeeze(sum(AllXYT,3));
All.TimeX = squeeze(sum(AllXYT,2))';
All.TimeY = squeeze(sum(AllXYT,1))';

All.X = ((1:xbins)-0.5) * areabin / ExpData.PixelSize;
All.Y = ((1:ybins)-0.5) * areabin / ExpData.PixelSize;
All.Time = ((1:tbins)-0.5) * timebin;

All.SingleWormSize = singleworm;
All.SizeHist = [wms; nw];
All.ClumpFrxn = nclump;

%--------------------
% print figure - All 
densfig = findobj(get(0,'Children'),'Tag','Worm Density');
if isempty(densfig) densfig = figure; set(densfig,'Tag','Worm Density'); end
figure(densfig); clf;

n = 3; for i=1:n; 
    subplot(5,n+1,8+i); 
    imagesc(All.X,All.Y,squeeze(sum(AllXYT(:,:,1+round(end*(i-1)/n):round(i*end/n)),3))'); 
    colormap([flipud(gray(64));zeros(64,3)]); axis equal tight; title(sprintf('%d/%d',i,n)); if i==1 ylabel('y(mm)'); end
end; subplot(5,n+1,8+(n+1)); imagesc(All.X,All.Y,All.XY'); axis equal tight; title('All'); xlabel('x(mm)');
hy = subplot(5,4,1:3); imagesc(All.Time,All.Y,All.TimeY'); colormap([flipud(gray(64));zeros(64,3)]); ylabel('y(mm)'); 
title(FileData(FileNum).SegFile,'Interpreter','none');
hx = subplot(5,4,5:7); imagesc(All.Time,All.X,All.TimeX'); colormap([flipud(gray(64));zeros(64,3)]); ylabel('x(mm)'); xlabel('Time (min)');
subplot(5,4,4); barh(All.Y,sum(All.XY,1)); ylim(get(hy,'YLim')); axis ij; title('(All animals)');
subplot(5,4,8); barh(All.X,sum(All.XY,2)); ylim(get(hx,'YLim')); axis ij;
warning('off','MATLAB:log:logOfZero');
warning('off','MATLAB:divideByZero');
subplot(5,4,13:14); bar(sizebins/singleworm,log10(hist(AllData(:,2),sizebins)));grid on; 
    title(sprintf('Object Size (1 worm = %0d pix^2)',round(singleworm))); 
    set(gca,'XTick',1/2:max(sizebins)/singleworm);
    ylabel('log10(frames)'); xlabel('Animals per object');
orient tall
saveas(densfig,strrepl(FileData(FileNum).SegFile,'_seg.mat','_density_1_all.pdf'));

%----------------------------------------------%
%               SEGMENT ANALYSIS               %
%----------------------------------------------%

disp([datestr(now),': Analyzing Turn Data']);
[TurnData,SegData,SegPos,SegDir,segidx] = getTurnData(Tracks);
disp([datestr(now),': Turn analysis complete']);

%----------------------------------------------%
%                FRAME ANALYSIS                %
%----------------------------------------------%

% Correlate original tracks with joined tracks  
OTrack = [Tracks.OriginalTrack];
NumAnimals = max(OTrack);

TrackData = [];
for Animal = 1:NumAnimals
    tracklist = find(OTrack == Animal);

    AnimalFr = [Tracks(tracklist).Frames]';
    AnimalCode = [Tracks(tracklist).Code]';
    AnimalRev = [Tracks(tracklist).Reverse]';
    AnimalDir = [Tracks(tracklist).PathAngle]';
    
    AnimalDirUD = -(AnimalDir > 30 & AnimalDir < 150) + (AnimalDir > 210 & AnimalDir < 330) + 2*isnan(AnimalDir);
    AnimalSpeed = [Tracks(tracklist).Distance]' / ExpData.PixelSize * ExpData.FrameRate; % mm/s
    AnimalPath = [];
    for tr = tracklist
        AnimalPath = [AnimalPath; Tracks(tr).Path];
    end
    
    Stall = (AnimalCode == 9);
    Omega = (AnimalCode == 7);
    NonOmega = (AnimalCode >=1 & AnimalCode <=4) + (AnimalCode == 6);
    AnimalBeh = 2*(NonOmega & ~AnimalRev) + 3*Stall + 4*(NonOmega & AnimalRev) ...
                + 5*(Omega & ~AnimalRev) + 6*(Omega & AnimalRev);
            
    AnimalSpeedFP = AnimalSpeed;
    AnimalSpeedFP(find(AnimalBeh ~= 2 & AnimalBeh ~= 3)) = NaN;  % speed for forward (2) and stall (3) only
    
    TrackData = [TrackData; ones(size(AnimalBeh))*Animal, AnimalFr/ExpData.FrameRate/60, AnimalPath, AnimalBeh, AnimalSpeedFP, AnimalDirUD, AnimalDir]; 
end

% Bin numbers (tbins, xbins,etc) should be the same as before...
TimeBin = ceil(TrackData(:,2)/timebin); tbins = max(TimeBin);
XBin = ceil(TrackData(:,3)/areabin); xbins = max(XBin);
YBin = ceil(TrackData(:,4)/areabin); ybins = max(YBin);
validxy = find(XBin > 0 & YBin > 0);
TrackData = TrackData(validxy,:);
XBin = XBin(validxy);
YBin = YBin(validxy);
clear validxy
DirY = TrackData(:,7)+2;  % -1,0,+1 --> 1,2,3

TrackXYTD = zeros(xbins, ybins, tbins, s, 'single');
tic;
for i=1:size(TrackData,1)
    TrackXYTD(XBin(i),YBin(i),TimeBin(i),DirY(i)) = TrackXYTD(XBin(i),YBin(i),TimeBin(i),DirY(i)) + 1;  % count number of animal-frames 
end

Track.XYTimeD = TrackXYTD;
t = toc; disp(['TrackXYTD analysis: ',num2str(t),' s']);

Track.XYTime = squeeze(sum(TrackXYTD,4));
Track.XY = squeeze(sum(Track.XYTime,3));
Track.TimeX = squeeze(sum(Track.XYTime,2))';
Track.TimeY = squeeze(sum(Track.XYTime,1))';
Track.YTimeD = squeeze(sum(TrackXYTD,1));

Track.X = ((1:xbins)-0.5) * areabin / ExpData.PixelSize;
Track.Y = ((1:ybins)-0.5) * areabin / ExpData.PixelSize;
Track.Time = ((1:tbins)-0.5) * timebin;

%--------------------
% print figure - Tracked 
DensityPlot(Track,FileData(FileNum).SegFile);
saveas(densfig,strrepl(FileData(FileNum).SegFile,'_seg.mat','_density_4_track.pdf'));

%----------------------------------
BehaviorList = {2, 3, 4, 5, 6, [5 6]}; % 0 is forward turn
BehaviorText = {'Forward Run','Pause','Short Reversal','Pirouette Reverse','Pirouette Forward','Pirouette FR'};

for b = 1:length(BehaviorList)

    Behavior = cell2mat(BehaviorList(b));
    BehXYTD = zeros(xbins, ybins, tbins, s, 'single');
    tic;
    for i=1:size(TrackData,1)
        BehXYTD(XBin(i),YBin(i),TimeBin(i),DirY(i)) = BehXYTD(XBin(i),YBin(i),TimeBin(i),DirY(i)) + sum(TrackData(i,5) == Behavior);  % count number of animal-frames with specific behavior
    end
    BehXYT = squeeze(sum(BehXYTD,4));
    
    Beh(b).XYTime = BehXYT;
    Beh(b).XY = squeeze(sum(BehXYT,3));
    Beh(b).TimeX = squeeze(sum(BehXYT,2))';
    Beh(b).TimeY = squeeze(sum(BehXYT,1))';
    Beh(b).Code = Behavior;
    
    Beh(b).XYTimeD = BehXYTD;
    Beh(b).XYD = squeeze(sum(BehXYTD,3));
    Beh(b).XTimeD = squeeze(sum(BehXYTD,2));
    Beh(b).YTimeD = squeeze(sum(BehXYTD,1));
    Beh(b).Code = Behavior;
    Beh(b).Text = BehaviorText(b);
    
    Beh(b).X = Track.X;
    Beh(b).Y = Track.Y;
    Beh(b).Time = Track.Time;
    t = toc; disp(['Behavior analysis ',num2str(b),' [',num2str(Behavior),']: ',num2str(t),' s']);

end

%--------------------
% print figure - behavior
for b = 1:length(BehaviorList)
    DensityPlot(Beh(b),FileData(FileNum).SegFile);
    saveas(densfig,strrepl(FileData(FileNum).SegFile,'_seg.mat',sprintf('_density_beh%d_%s.pdf',b,char(Beh(b).Text))));
end

%-----------------------------------
SpeedXYTD = zeros(xbins, ybins, tbins, s, 'single');
tic;
validindex = find(TrackData(:,5) == 2 | TrackData(:,5) == 3); % find forward (2) and pause (3) states

for i=1:length(validindex)
    SpeedXYTD(XBin(validindex(i)),YBin(validindex(i)),TimeBin(validindex(i)),DirY(validindex(i))) = ...
        SpeedXYTD(XBin(validindex(i)),YBin(validindex(i)),TimeBin(validindex(i)),DirY(validindex(i))) + ...
        TrackData(validindex(i),6);  % sum speed of forward+pause
end

Speed.wt.XYTimeD = (Beh(1).XYTimeD + Beh(2).XYTimeD) > 0;
Speed.FP.XYTimeD = SpeedXYTD ./ ((Beh(1).XYTimeD + Beh(2).XYTimeD)+ ~Speed.wt.XYTimeD);
t = toc; disp(['SpeedXYTD analysis: ',num2str(t),' s']);

Speed.wt.XYTime = sum(Speed.wt.XYTimeD,4) > 0;
Speed.FP.XYTime = sum(Speed.FP.XYTimeD,4) ./ (sum(Speed.wt.XYTimeD,4) + (Speed.wt.XYTime == 0));

Speed.wt.XY = squeeze(sum(Speed.wt.XYTime,3));
Speed.FP.XY = squeeze(sum(Speed.FP.XYTime,3)) ./ (Speed.wt.XY + (Speed.wt.XY == 0));
Speed.wt.TimeX = squeeze(sum(Speed.wt.XYTime,2))';
Speed.FP.TimeX = squeeze(sum(Speed.FP.XYTime,2))' ./ (Speed.wt.TimeX + (Speed.wt.TimeX == 0));
Speed.wt.TimeY = squeeze(sum(Speed.wt.XYTime,1))';
Speed.FP.TimeY = squeeze(sum(Speed.FP.XYTime,1))' ./ (Speed.wt.TimeY + (Speed.wt.TimeY == 0));

Speed.wt.XYD = squeeze(sum(Speed.wt.XYTimeD,3));
Speed.FP.XYD = squeeze(sum(Speed.FP.XYTimeD,3)) ./ (Speed.wt.XYD + (Speed.wt.XYD == 0));
Speed.wt.XTimeD = squeeze(sum(Speed.wt.XYTimeD,2));
Speed.FP.XTimeD = squeeze(sum(Speed.FP.XYTimeD,2)) ./ (Speed.wt.XTimeD + (Speed.wt.XTimeD == 0));
Speed.wt.YTimeD = squeeze(sum(Speed.wt.XYTimeD,1));
Speed.FP.YTimeD = squeeze(sum(Speed.FP.XYTimeD,1)) ./ (Speed.wt.YTimeD + (Speed.wt.YTimeD == 0));

Speed.FP.Time = Track.Time;
Speed.FP.X = Track.X;
Speed.FP.Y = Track.Y;

%--------------------
% print figure - speed

densfig = findobj(get(0,'Children'),'Tag','Worm Density');
if isempty(densfig) densfig = figure; set(densfig,'Tag','Worm Density'); end
figure(densfig); clf;
dirlabel = {'up','horiz','down'};

hy = subplot(5,4,1:3); imagesc(Speed.FP.Time,Speed.FP.Y,Speed.FP.TimeY'); colormap([1 1 1; jet(64)]); ylabel('y(mm)'); 
title(FileData(FileNum).SegFile,'Interpreter','none');
hx = subplot(5,4,5:7); imagesc(Speed.FP.Time,Speed.FP.X,Speed.FP.TimeX'); colormap([1 1 1; jet(64)]); ylabel('x(mm)'); xlabel('Time (min)');
subplot(5,4,4); barh(Speed.FP.Y,sum(Speed.FP.XY .* Speed.wt.XY,1) ./ sum(Speed.wt.XY,1)); ylim(get(hy,'YLim')); axis ij; title('Fwd Speed (mm/s)'); 
subplot(5,4,8); barh(Speed.FP.X,sum(Speed.FP.XY .* Speed.wt.XY,2) ./ sum(Speed.wt.XY,2)); ylim(get(hx,'YLim')); axis ij; 
n = 3; for j=1:3; for i=1:n; 
    subplot(5,n+1,8+((j-1)*(n+1))+i); 
    imagesc(Speed.FP.X,Speed.FP.Y,squeeze(msum(Speed.FP.XYTimeD(:,:,1+round(end*(i-1)/n):round(i*end/n),j) .* ...
                                               Speed.wt.XYTimeD(:,:,1+round(end*(i-1)/n):round(i*end/n),j),3:4) ./ ...
                                          msum(Speed.wt.XYTimeD(:,:,1+round(end*(i-1)/n):round(i*end/n),j),3:4))'); 
    colormap([1 1 1; jet(64)]); axis equal tight; title(sprintf('%d/%d %s',i,n,char(dirlabel(j)))); if i==1 ylabel('y(mm)'); end; if j==3 xlabel('x(mm)'); end;
    end;
    subplot(5,n+1,8+j*(n+1)); 
    imagesc(Speed.FP.X,Speed.FP.Y,squeeze(msum(Speed.FP.XYTimeD(:,:,:,j) .* ...
                                               Speed.wt.XYTimeD(:,:,:,j),3:4) ./ ...
                                          msum(Speed.wt.XYTimeD(:,:,:,j),3:4))'); 
    axis equal tight; title(sprintf('All %s',char(dirlabel(j)))); if j==3 xlabel('x(mm)'); end;
end 

orient tall
saveas(densfig,strrepl(FileData(FileNum).SegFile,'_seg.mat','_density_3_speed.pdf'));

%------------------------
% print figure - summary

figure(densfig); clf;
hy = subplot(7,4,1:3); imagesc(Track.Time,Track.Y,Track.TimeY'); colormap([flipud(gray(64));zeros(64,3)]); ylabel('y(mm)'); 
title(FileData(FileNum).SegFile,'Interpreter','none');
subplot(7,4,4); barh(Track.Y,sum(Track.XY,1)); ylim(get(hy,'YLim')); axis ij; 

linecolor = [1 0 0; .6 1 .6; 0 0 1];
numybins = sum(sum(Track.XY,2)>0);
n = 3; for j=[2 1 3]; for i=1:n+1;
    if i<=n timerange = 1+round(size(Track.XYTimeD,3)*(i-1)/n):round(i*size(Track.XYTimeD,3)/n); titlestr = sprintf('%d/%d',i,n);
    else timerange = 1:size(Track.XYTimeD,3); titlestr = 'All'; end
        
    subplot(7,n+1,4+i); hold on; title(titlestr);
        stairs(Track.Y,mean(msum(Track.XYTimeD(:,:,timerange,j),[1 4]) ./ repmat(msum(Track.XYTimeD(:,:,timerange,j),[1 2 4]),1,size(Track.XYTimeD,2))',2) * numybins, ...
        'Color',linecolor(j,:)); 
        axis xy; xlim(get(hy,'YLim')); ylim([0 4]); if i==1 ylabel('Residence'); end
    
    subplot(7,n+1,4+n+1+i); hold on;
        stairs(Speed.FP.Y,squeeze(msum(Speed.FP.XYTimeD(:,:,timerange,j) .* ...
                            Speed.wt.XYTimeD(:,:,timerange,j),[1 3 4]) ./ ...
                       msum(Speed.wt.XYTimeD(:,:,timerange,j),[1 3 4])),'Color',linecolor(j,:)); 
        axis xy; xlim(get(hy,'YLim')); ylim([0 0.4]); if i==1 ylabel('Fwd Speed (mm/s)'); end
        
    subplot(7,n+1,4+2*(n+1)+i); hold on;
        stairs(Track.Y,msum(Beh(3).XYTimeD(:,:,timerange,j),[1 3 4]) ./ msum(Track.XYTimeD(:,:,timerange,j),[1 3 4]), 'Color',linecolor(j,:)); 
        axis xy; xlim(get(hy,'YLim')); ylim([0 0.4]); if i==1 ylabel('Rev Prob.'); end

    if j~=2
        subplot(7,n+1,4+3*(n+1)+i); hold on;
            groupidx = findpattern(TurnData(:,[1,4:6]),[min(timerange)*60 max(timerange)*60; 4 4; 2 2; j j]'); % surf turns
            stairs(Track.Y,hist(TurnData(groupidx,3)/10,1:length(Track.Y)),'Color',linecolor(j,:)); 
            axis xy; xlim(get(hy,'YLim')); if i==1 ylabel('Surf Turns'); end

        subplot(7,n+1,4+4*(n+1)+i); hold on; exitdir = [2 3; 1 3; 1 2];
            groupidx = findpattern(TurnData(:,[1,4:6]),[min(timerange)*60 max(timerange)*60; 1 6; j j; exitdir(j,:)]'); % up/down turns
            stairs(Track.Y,hist(TurnData(groupidx,3)/10,1:length(Track.Y)),'Color',linecolor(j,:)); 
            axis xy; xlim(get(hy,'YLim')); if i==1 ylabel('Turns up/down'); end
    end
    
    subplot(7,n+1,4+5*(n+1)+i); hold on;
        groupidx = findpattern(TurnData(:,[1,4:5]),[min(timerange)*60 max(timerange)*60; 7 7; j j]'); % pirouette
        stairs(Track.Y,hist(TurnData(groupidx,3)/10,1:length(Track.Y)),'Color',linecolor(j,:)); 
        axis xy; xlim(get(hy,'YLim')); if i==1 ylabel('Pirouette'); end

end; end 

% add legends
legendh = axes('Position',[0.9 0.1 0.05 0.8],'box','off','visible','off');
texth(1) = text(0.5,0.6,'moving up (-y)','FontSize',8,'Color',[1 0 0],'Rotation',90);
texth(2) = text(0.8,0.6,'moving horiz','FontSize',8,'Color',[.6 1 .6],'Rotation',90);
texth(3) = text(1.1,0.6,'moving down (+y)','FontSize',8,'Color',[0 0 1],'Rotation',90);
texth(4) = text(0.5,0.3,'curve horiz --> up','FontSize',8,'Color',[1 0 0],'Rotation',90);
texth(5) = text(1.1,0.3,'curve horiz --> down','FontSize',8,'Color',[0 0 1],'Rotation',90);
texth(6) = text(0.5,0.1,'heading up','FontSize',8,'Color',[1 0 0],'Rotation',90);
texth(7) = text(0.8,0.1,'heading horiz','FontSize',8,'Color',[.6 1 .6],'Rotation',90);
texth(8) = text(1.1,0.1,'heading down','FontSize',8,'Color',[0 0 1],'Rotation',90);


orient tall
saveas(densfig,strrepl(FileData(FileNum).SegFile,'_seg.mat','_density_2_summary.pdf'));

PDFcomplete = false;
while ~PDFcomplete
    pause(1);
    PDFcomplete = exist(strrepl(FileData(FileNum).SegFile,'_seg.mat','_density_2_summary.pdf')) == 2;
end
    
%-----------
% Save Data 
%-----------

if SaveData
    FileInfo = FileData(FileNum);
    
    saveext = '_density';
    FileInfo.SaveFile = fullfile(FileInfo.PathName,[FileInfo.BaseName,saveext,Ext]);
    FileInfo.FileNum = [FileNum, NumFiles];
    
    ExpData.TimeBin = timebin;
    ExpData.AreaBin = areabin;

    disp([datestr(now),'-- Saving Binned Analysis File ',num2str(FileNum)]);
    save(FileInfo.SaveFile,'FileInfo','ExpData','TrackData','All','Track','Beh','Speed', ...
        'TurnData','SegData','SegPos','SegDir','segidx');
    disp([datestr(now),'--     ...saved.']);
    
    BinFileList = cat(1,BinFileList,{FileInfo.SaveFile});
end

disp('Assembling output pdfs...');
pause(10);
currentpath = pwd;
cd(FileData(FileNum).PathName);
onepdf(['*',FileData(FileNum).BaseName,'*.pdf'],[FileData(FileNum).BaseName,'_analysis.pdf']); 
cd(currentpath);
% onepdf(['"',FileData(FileNum).PathName,'\*',FileData(FileNum).BaseName,'*.pdf"'],['"',FileData(FileNum).PathName,'\',FileData(FileNum).BaseName,'_analysis.pdf"']);

end
%---------------------------------------------%
%                END BIG LOOP                 %
%---------------------------------------------%