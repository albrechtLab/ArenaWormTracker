%
% Track worms from video file (AVI, uncompressed or Indeo5 compressed). 
% Requires user input:
%   (1) Tracking options: range of video frames to track, fragments to
%   divide up the movie, video framerate, and arena size.
%   (2) Arena bounds: Draw a box.
%   (3) Vertical arena dimensions. Click upper and lower positions that
%   correspond to arena size.
% Saves track data to a file.
%
% USAGE:
%   File = ArenaTracker(moviefile,blankrect)
%
%   moviefile: input video filename
%   blankrect: position of rectangle to ignore
%   File: structure containing file info

%---------------------------- 
% Dirk Albrecht 
% Version 1.2 
% 29-Jan-2021
%---------------------------- 

function File = ArenaTracker(moviefile,blankrect);

%keep3('command','blankrect');       % clear variables
close(findobj('type','figure'));    % close all figures

%----------------------------------------
% Define Constants & Intialize variables
%----------------------------------------

arena = 16.1;                 % mm arena height (perpendicular to flow), typically 22 or 16.1 mm
FrameRate = 2;              % Framerate (Hz)
AutoSelectDyeTrack = 1;
TrackDye = 1;
bs = 5;                     % pixels for background and dye levels
DefaultFragmentLength = 1800; % frames

%-----------------
% Set screen info
%-----------------
WTFigH = findobj('Tag', 'WTFIG');
if isempty(WTFigH)
    WTFigH = figure('NumberTitle', 'off', 'Tag', 'WTFIG');
else
    figure(WTFigH);
end
scrsz = get(0,'ScreenSize');
  
%----------------
% Get File info 
%----------------
MovieNum = 1;
UserPrompt = ['Select Video File For Analysis:'];
[FileName, PathName] = uigetfile({'*.avi;*.mp4', 'Video files (*.avi, *.mp4)'}, UserPrompt);
cd(PathName);
Name(MovieNum).MovieFile = fullfile(PathName,FileName);
[~,~,ext] = fileparts(FileName);
Name(MovieNum).TrackFile = strrepl(Name(MovieNum).MovieFile, ext,'.mat');

%-----------------
% Get User Input 
%-----------------
MovieName = Name(MovieNum).MovieFile;
[p,ShortMovieName] = fileparts(MovieName);

TrackSettingsName = fullfile(p,[ShortMovieName,'_TrackSettings.mat']);
CollectTrackingSettings = true;
if exist(TrackSettingsName) == 2
    ButtonName = questdlg('Use existing tracking settings?','','Yes', 'No', 'Yes');
    if strcmp(ButtonName,'Yes')
        try
            load(TrackSettingsName);
            File = FileSettings;

            % check filenames
            [pn1,fn1] = fileparts(Name(MovieNum).MovieFile); % chosen
            [pn2,fn2] = fileparts(File(MovieNum).MovieFile); % from settings file
            if strcmp(fn1,fn2) % if filenames match
                if ~strcmp(pn1,pn2) % if pathnames don't match
                    File(MovieNum).MovieFile = Name(MovieNum).MovieFile;
                    File(MovieNum).TrackFile = Name(MovieNum).TrackFile;
                end
                CollectTrackingSettings = false;
            else
                disp('Filenames don''t match... reselect manually.');
            end

        catch
            disp('Error loading settings.  Select manually...');
        end
    end
end

if CollectTrackingSettings

    % show first frame
    MovieObj = VideoReader(MovieName);
    Mov = read(MovieObj, 1);
    Mov = Mov(:,:,3);
    MovSize = size(Mov);
    figure(WTFigH); clf; imshow(imadjust(Mov,[],[],2));
    set(gcf,'MenuBar','none');
    set(gcf,'Toolbar','none');
    set(gcf,'Position',[150 100 1020 900]);

    %--------get information----------
    buttony = 10;

    FrameNum = MovieObj.NumberOfFrames;
    FrameStart = 1;
    NumFragments = max(round(FrameNum/DefaultFragmentLength),1);

    uicontrol('Style','text','Position',[230 buttony 40 20],'String','Frames');
    h2a = uicontrol('style','edit','Position',[270 buttony 50 20],'String',num2str(FrameStart));
    h2 = uicontrol('style','edit','Position',[320 buttony 50 20],'String',num2str(FrameNum));
    h2b = uicontrol('style','edit','Position',[380 buttony 20 20],'String',num2str(NumFragments));

    h3 = uicontrol('style','checkbox','Position',[420 buttony 60 20],'String','Track Dye','Value',TrackDye);

    uicontrol('Style','text','Position',[500 buttony 60 20],'String','Arena(mm)');
    h4 = uicontrol('style','edit','Position',[560 buttony 40 20],'String',num2str(arena));

    uicontrol('Style','text','Position',[640 buttony 60 20],'String','FrameRate');
    h5 = uicontrol('style','edit','Position',[700 buttony 30 20],'String',num2str(FrameRate));

    %--------------------------

    set(WTFigH,'Name',['Movie ',num2str(MovieNum),': ',ShortMovieName]);
    txt = 'Click background (PDMS)'; title(txt);
    label = text(MovSize(2)/2,MovSize(1)*0.4,txt,'FontSize',18,'HorizontalAlignment','center');
    [X,Y] = ginput(1); File(MovieNum).BgBox = [X-bs Y-bs 2*bs 2*bs]; 
    rectangle('Position',File(MovieNum).BgBox); text(X+10,Y,'bg');

    ar = 1; a = 1;
    % get tracking regions
    txt = 'SELECT TRACKING AREA: click inside to confirm, outside to redo.';
    title(txt); set(label,'String',txt);
    success = 0; while success == 0
        box = getrect(gca); box = box+(box==0); h = rectangle('Position',box); set(h,'EdgeColor',[1,0,0]);
        [X,Y,button] = ginput(1); 
        if button == 1 & (X-box(1) >= 0 & X-box(1) <= box(3) & Y-box(2) >= 0 & Y-box(2) <= box(4))
            File(MovieNum).Arena(ar).TrackBox = box+(box==0);
            set(h,'EdgeColor',[0,0,1]); 
            success = true;
        end
    end

    File(MovieNum).NumArenas = ar;

    TrackDye = get(h3,'Value');
    % get dye information
    if TrackDye
        if AutoSelectDyeTrack
            box = File(MovieNum).Arena(a).TrackBox;
            File(MovieNum).Arena(a).UpBox = [box(1)-15-bs box(2)+box(4)/2-2*bs 2*bs 4*bs]; 
                rectangle('Position',File(MovieNum).Arena(a).UpBox); text(X+10,Y,[num2str(a),'up']);
            File(MovieNum).Arena(a).DnBox = [box(1)+box(3)+15-bs box(2)+box(4)/2-2*bs 2*bs 4*bs]; 
                rectangle('Position',File(MovieNum).Arena(a).DnBox); text(X+10,Y,[num2str(a),'dn']);
        else
            txt = ['Pick upstream and downstream points for Arena: ',num2str(a)];
            title(txt); set(label,'String',txt);
            [X,Y] = ginput(1); File(MovieNum).Arena(a).UpBox = [X-bs Y-2*bs 2*bs 4*bs]; 
                rectangle('Position',File(MovieNum).Arena(a).UpBox); text(X+10,Y,[num2str(a),'up']);
            [X,Y] = ginput(1); File(MovieNum).Arena(a).DnBox = [X-bs Y-2*bs 2*bs 4*bs]; 
                rectangle('Position',File(MovieNum).Arena(a).DnBox); text(X+10,Y,[num2str(a),'dn']);
        end
    else
        File(MovieNum).Arena(a).UpBox = [1 1 1 1];
        File(MovieNum).Arena(a).DownBox = [1 1 1 1];
    end

    % get scale info
    txt = ['GET SCALING: pick points on top and bottom arena edges = ',get(h4,'String'),'mm'];
    title(txt); set(label,'String',txt);
    [X,Y] = ginput(2); 
    arenapix = abs(Y(2)-Y(1)); 

    % extra rectangular region for exclusion 
    if exist('blankrect')
        txt = 'SELECT EXTRA EXCLUSION REGION';
        title(txt); set(label,'String',txt);
        box = getrect(1); box = box+(box==0); 
        blankrect = round(max(box,1));
    end

    File(MovieNum).MovieFile = Name(MovieNum).MovieFile;
    File(MovieNum).TrackFile = Name(MovieNum).TrackFile;

    File(MovieNum).TrackFrames = str2num(get(h2,'String'));
    File(MovieNum).StartFrame = str2num(get(h2a,'String'));
    File(MovieNum).DyeTracked = get(h3,'Value');
    File(MovieNum).ArenaSize = str2num(get(h4,'String'));
    File(MovieNum).PixelSize = arenapix / File(MovieNum).ArenaSize;
    File(MovieNum).ImageSize = MovSize;
    File(MovieNum).FrameRate = str2num(get(h5,'String'));
    File(MovieNum).Fragments = str2num(get(h2b,'String'));
    File(MovieNum).StartFragment = 1;
    File(MovieNum).CompletedFragments = zeros(1,File(MovieNum).Fragments);
    File(MovieNum).ClaimedFragments = zeros(1,File(MovieNum).Fragments);
    File(MovieNum).Stitched = 0;

    hi = File(MovieNum).TrackFrames;
    lo = File(MovieNum).StartFrame;
    NumFragments = File(MovieNum).Fragments;
    perfrag = round((hi-lo)/NumFragments/60)*60;
    frfr = [lo+(0:(NumFragments-1))*perfrag]';
    frfr = [frfr, [frfr(2:NumFragments)-1; hi]];

    File(MovieNum).FragmentFrames = frfr;
    
    disp([datestr(now),' Image Data Collected ']);
    FileSettings = File(MovieNum);
    save(TrackSettingsName,'FileSettings');
end

%---------------
% Start Tracker
%---------------

FrameNum = File(MovieNum).TrackFrames;
Start = File(MovieNum).StartFrame;
ImageSize = File(MovieNum).ImageSize;
FrameRate = File(MovieNum).FrameRate;
PixelSize = File(MovieNum).PixelSize;
FragmentFrames = File(MovieNum).FragmentFrames;
NumArenas = File(MovieNum).NumArenas;
TrackDye = File(MovieNum).DyeTracked;
NumFragments = File(MovieNum).Fragments;
       
TrackName = File(MovieNum).TrackFile;
MovieName = File(MovieNum).MovieFile;
    
[pathname,filename,ext] = fileparts(TrackName);
partsfolder = '\trackfile parts\';
FragmentSaveNames = {}; 
for i = 1:File(MovieNum).Fragments;
    FragmentSaveNames = [FragmentSaveNames; {[pathname,partsfolder,filename,'_part',num2str(i),ext]}]; 
end
if exist([pathname,partsfolder]) == 0 mkdir([pathname,partsfolder]); end;

[p,ShortMovieName] = fileparts(MovieName);

Tracks = [];
AllData = [];
DyeData = [];
DyeData.Time = [];
DyeData.Background = [];
ExpData = [];

%Video Mask:  remove bright or dim pixels (deals with timer)
MovieObj = VideoReader(MovieName);
Mov = read(MovieObj, MovieObj.NumberOfFrames);
MovData = Mov(:,:,3);
MaskPix = 10; MaskBorder = 30;
Mask = (MovData <= MaskPix | MovData >= 255-MaskPix);
Mask = imdilate(Mask,strel('square',3));
Mask(MaskBorder:ImageSize(1)-MaskBorder, MaskBorder:ImageSize(2)-MaskBorder) = 0;

if exist('blankrect') 
    Mask(blankrect(2)+(0:blankrect(4)), blankrect(1)+(0:blankrect(3))) = 1;
end

tic;

%----------------------
% ANALYZE FRAGMENTS
%----------------------
    
for Fragment = 1:File(MovieNum).Fragments
   FragmentTracker(File,FragmentSaveNames,FragmentFrames,Fragment,Mask,TrackDye,WTFigH); 
end
    
% Save Combined File
StitchTracks(FragmentSaveNames,TrackName);

% Copy files back to original location
disp([datestr(now),': Copying tracks file from local directory']);
[s,mess,messid] = copyfile(TrackName,File(MovieNum).TrackFile);
disp([datestr(now),': Complete... ',mess]);

disp([datestr(now),': Copying avi.mat file(s) from local directory']);
[s,mess,messid] = copyfile([MovieName,'*.mat'],fileparts(File(MovieNum).TrackFile));
disp([datestr(now),': Complete... ',mess]);

