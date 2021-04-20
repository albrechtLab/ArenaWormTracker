%
% Identifies instantaneous behavioral state from worm tracks and
% morphological data from ArenaTracker.
% Saves behavioral data to a file (*_seg.mat)
%
% USAGE:
%   SaveList = SegmentTracks(FileName)
%
%   FileName: single filename of tracked data or cell array of multiple
%               filenames. Select with user input if none given.
%   SaveList: cell array of output filename(s).

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 30-Mar-2011 14:17:48 
%---------------------------- 

function SaveList = SegmentTracks(FileName)
    
if nargin < 1
    FileName = {};
end

SaveList = {};

% Initialize Segment Analysis Settings
SegmentationSettings;

if length(FileName) == 0
    % Get track data for analysis
    % --------------------------
    [FileName, PathName] = uigetfile('*.mat','Select Track File(s) For Analysis','MultiSelect','on');
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
disp(['Batch Segmenting ',num2str(NumFiles),' tracks files.']);

for fnum = 1:NumFiles
%--------------------------------------------------------------------------
% BIG LOOP FOR EACH FILE 
%--------------------------------------------------------------------------

if exist('PathName')
    FullName = fullfile(PathName, char(FileName(fnum)));
else
    FullName = char(FileName(fnum));
end

disp(sprintf('File %d of %d: %s',fnum,NumFiles,FullName));

[pathname,filename,ext] = fileparts(FullName);  SaveName = [filename,'_seg',ext];
FullSaveName = fullfile(pathname, SaveName);
BehSaveName = [FullSaveName(1:end-8),'_beh.mat'];
 
%------------
% Load Data
%------------
clear('Tracks','background','DyeData','ExpData','File','AllData');

folderdate = []; founddate = strfind(FullName,'200'); if length(founddate)>0 folderdate = FullName(founddate+(0:7)); end
disp([datestr(now),': Loading File #',num2str(fnum),': [',folderdate,'] ',filename]);

load(FullName,'Tracks','background','DyeData','ExpData','File','AllData');

if ~exist('Tracks')
    disp('no Tracks variable found. skipping...');
    continue
end
disp(['...loaded ',datestr(now)]);

AllData = single(AllData);

if (exist('ExpData') && length(ExpData) > 0)
    framerate = ExpData.FrameRate;
    pixelsize = ExpData.PixelSize;
else
    framerate = 2;
    pixelsize = 1;
    disp('WARNING: NO SCALING DATA');
end
Settings.FrameRate = framerate;
Settings.PixelSize = pixelsize;

if isfield(Tracks,'Segment')
    ButtonName=questdlg('Segments have already been analyzed.','','Reanalyze', 'Stop', 'Reanalyze');
    switch ButtonName
        case 'Stop'
            disp('ending.'),
            return;
    end
end

TrackArena = ones(1, length(Tracks));
NumArenas = 1;
%---------------------------
numtracks = length(Tracks);
trackstats = []; 
totalfr = 0;
sprg = -1;
tic;

for tr = 1:length(Tracks)
%--------------------------------------------------------------------------
%  LOOP FOR EACH TRACK
%--------------------------------------------------------------------------

    [Segment, SegAnalysis] = SegmentTrack(Tracks(tr),Settings);
    trackbox = DyeData.Arena(TrackArena(tr)).TrackBox;

    TrackAnalysis = AnalyzeTrack(Segment,SegAnalysis,trackbox);

    reversal = (TrackAnalysis.FRP == 2);

    %------------------------------
    % reorganize data for plotting
    %------------------------------
    % segment codes:
    %   1 - fwd
    %   2 - fwd w/ revs
    %   3 - lane chg
    %   4 - small curve (60deg)
    %   5 - short reversal
    %   6 - long curve / loopy movement
    %   7 - pirouette (sharp curve following reversal)
    %   8 - unknown/erratic movement
    %   9 - pause
    %  10 - out of bounds / collision

    distance = []; pathang = []; pathangvel = [];
    for seg = 1:length(Segment)
        distance(Segment(seg).RealFrames) = Segment(seg).Distance';
        pathang(Segment(seg).RealFrames) = Segment(seg).PathAngle';
        pathangvel(Segment(seg).RealFrames) = Segment(seg).PathAngVel';
    end

    distance(SegAnalysis.StallFr) = 0;
    pathang(SegAnalysis.StallFr) = NaN;
    pathangvel(SegAnalysis.StallFr) = NaN;

    spd = distance / ExpData.PixelSize * ExpData.FrameRate .* (1-2*reversal);

    t = toc;
    trlength = length(TrackAnalysis.Code);
    numsegments = length(Segment);
    totalfr = totalfr + trlength;

    status = sprintf('Track: %d/%d [%d fr] %d @ %d fps]',tr,numtracks,trlength,totalfr,round(totalfr/t));
    sprg = showprog(status,sprg);

    trackstats = [trackstats; trlength, numsegments];

    Tracks(tr).Code = TrackAnalysis.Code;
    Tracks(tr).Distance = distance;
    Tracks(tr).Speed = spd;
    Tracks(tr).PathAngle = pathang;
    Tracks(tr).PathAngVel = pathangvel;
    Tracks(tr).Reverse = reversal;
    Tracks(tr).OriginalDistance = SegAnalysis.OriginalDistance;
    Tracks(tr).Segment = Segment;
    Tracks(tr).Stall = SegAnalysis.StallFr;
    Tracks(tr).NoStall = SegAnalysis.NoStallFr;
    Tracks(tr).Omega = SegAnalysis.Omega;
    Tracks(tr).FTurnCCW = SegAnalysis.FTurnCCW;
    if size(Tracks(tr).FTurnCCW,1) > 0 Tracks(tr).FTurnCCW = Tracks(tr).FTurnCCW(find(TrackAnalysis.Code(Tracks(tr).FTurnCCW(:,1)) ~= 10),:); end
    Tracks(tr).FTurnCW = SegAnalysis.FTurnCW;
    if size(Tracks(tr).FTurnCW,1) > 0 Tracks(tr).FTurnCW = Tracks(tr).FTurnCW(find(TrackAnalysis.Code(Tracks(tr).FTurnCW(:,1)) ~= 10),:); end

    Tracks(tr).TurnOmega = TrackAnalysis.TurnOmega;
    Tracks(tr).Beh = TrackAnalysis.Beh;
end
         
ExpData.NumArenas = NumArenas;
if ~isfield(ExpData,'FrameRate') ExpData.FrameRate = 2; end

Arena = 1;
ExpData.Arena = Arena;

    %-----------------------
    % Link Tracks together
    %-----------------------
    linkoutput = LinkTracks(Tracks,0.5*ExpData.PixelSize/ExpData.FrameRate,Inf);  %Max avg velocity = 0.5mm/s
    [Tracks.OriginalTrack] = deal(linkoutput.OriginalTrack);
    ExpData.Animals = max(struct2mat(1,linkoutput,[],{'OriginalTrack'}));

    SegSettings = Settings;
    
    save(FullSaveName,'Tracks','background','DyeData','ExpData','File','AllData','trackstats','SegSettings');

    %Save condensed data
    disp('Condensing behavior data...'); 
    
    BehSaveName = [FullSaveName(1:end-8),'_beh.mat'];
    tfields = fieldnames(Tracks); 
    rfields = tfields(find(~strcmp(tfields,'Beh') & ...
                           ~strcmp(tfields,'Frames') & ...
                           ~strcmp(tfields,'OriginalTrack') & ...
                           ~strcmp(tfields,'PathAngle') & ...
                           ~strcmp(tfields,'Speed')));
    for tr=1:length(Tracks)
        Tracks(tr).X = Tracks(tr).Path(:,1)'; 
        Tracks(tr).Y = Tracks(tr).Path(:,2)'; 
    end                   
    Tracks = rmfield(Tracks,rfields);
    Tracks = singleStruct(Tracks);

    save(BehSaveName,'Tracks');    
    disp(['done. ',datestr(now)]);
    
    SaveList = cat(1,SaveList,{FullSaveName});
end

%---------------------end big loop for each file
end

