%
% Track worms for a movie fragment.
%
% USAGE:
%   FragmentTracker(File,FragmentSaveNames,FragmentFrames,Fragment,Mask,TrackDye,WTFigH)
%
%   File: file structure from ArenaTracker
%   FragmentSaveNames: list of save names
%   FragmentFrames: list of frames for all fragments
%   Fragment: current fragment number
%   Mask: video mask pattern
%   TrackDye: flag to track dye concentration
%   WTFigH: figure handle

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 24-Mar-2011 22:01:19 
%---------------------------- 
function FragmentTracker(File,FragmentSaveNames,FragmentFrames,Fragment,Mask,TrackDye,WTFigH)

MovieNum = 1;
NumArenas = 1;
MinWormArea = 0.5;          % Min area for object to be a valid worm
MaxWormArea = 1.7;          % Max area for single worm
MinDistance = 15;           % Min Distance for connecting a new worm to an existing track (in pixels)
SizeChangeThreshold = 100;  % Max size change between frames (in pixels)
MinTrackLength = 10;        % Min Length of valid track (in frames)    
PlotFrameRate = 60;         % Display tracking results every 'PlotFrameRate' frames

% Get background
background = getbackground(File(MovieNum).MovieFile,FragmentFrames(Fragment,1),FragmentFrames(Fragment,2),50);
disp([datestr(now),' Background Processed for Movie ',num2str(MovieNum),', fragment ',num2str(Fragment)]);
fprintf('\n');

% Initialize variables
Tracks = [];
AllData = [];
DyeData.Time = []; 
DyeData.Background = [];
for a = 1:NumArenas
    DyeData.Arena(NumArenas).Up = []; DyeData.Arena(NumArenas).Down = [];
end

Level = 0;
Settings = [];

MovieObj = VideoReader(File(MovieNum).MovieFile);

for Frame = FragmentFrames(Fragment,1):FragmentFrames(Fragment,2)
   
    % Get Frame
    Mov = read(MovieObj, Frame);

    %subtract the background from the frame
    Movdiv = 1 - min(double(Mov(:,:,1)+1)./(double(background)+1),1);
    %Movdiv(TimerY,TimerX) = 0;  % clear timer region
    Movdiv = Movdiv .* ~Mask;

    % calculate threshold and animal size
    if Level == 0
        [Level,AnimalPix] = AutoThreshold(Movdiv); 
        Settings = [Settings; Frame, Level, AnimalPix];
        disp(sprintf('Auto-threshold: %0.3f (Animal = %d pix; %0.3f mm^2)',Level,AnimalPix,AnimalPix / File.PixelSize^2));
        OLevel = Level;
    end

    % Filter and threshold video frame
    F2 = VideoFilter(Movdiv);
    BW = im2bw(F2,Level);

    % Segment objects    
    [L,NUM] = bwlabel(BW);
    if NUM > 300 
        % Reset autothreshold if too many objects
        [Level,AnimalPix] = AutoThreshold(Movdiv);
        Settings = [Settings; Frame, Level, AnimalPix]; disp(['(',num2str(Level),')']);
        BW = im2bw(F2,Level);
        [L,NUM] = bwlabel(BW);
    end

    STATS = regionprops(L, {'Area', 'Centroid', 'Eccentricity', 'MajorAxisLength', ...
        'MinorAxisLength', 'Orientation', 'Image', 'BoundingBox'});

    % Identify all worms by size, get their centroid coordinates
    WormIndices = find([STATS.Area] > MinWormArea*AnimalPix & [STATS.Area] < MaxWormArea*AnimalPix);

    NumWorms = length(WormIndices);
    if isempty(NumWorms)
        break
    end

    WormCentroids = [STATS(WormIndices).Centroid];
    WormCoordinates = [WormCentroids(1:2:2*NumWorms)', WormCentroids(2:2:2*NumWorms)'];
    WormSizes = [STATS(WormIndices).Area];
    WormEccentricities = [STATS(WormIndices).Eccentricity];
    WormMajorAxes = [STATS(WormIndices).MajorAxisLength];
    WormMinorAxes = [STATS(WormIndices).MinorAxisLength];
    WormOrientation = [STATS(WormIndices).Orientation];
    WormBox = [STATS(WormIndices).BoundingBox];
    WormBoundingBox = [WormBox(1:4:4*NumWorms)', WormBox(2:4:4*NumWorms)'];
    WormImage = [];
    for wi = 1:length(WormIndices)
        WormImage(wi).Image = STATS(WormIndices(wi)).Image;
    end

    for i=1:length(STATS) 
        AllData = [AllData; Frame, STATS(i).Area, STATS(i).Centroid, STATS(i).Eccentricity];
    end

    % Track worms
    % -----------
    if ~isempty(Tracks)
        ActiveTracks = find([Tracks.Active]);
    else
        ActiveTracks = [];
    end

    % Update active tracks with new coordinates
    for i = 1:length(ActiveTracks)
        DistanceX = WormCoordinates(:,1) - Tracks(ActiveTracks(i)).LastCoordinates(1);
        DistanceY = WormCoordinates(:,2) - Tracks(ActiveTracks(i)).LastCoordinates(2);
        Distance = sqrt(DistanceX.^2 + DistanceY.^2);
        [MinVal, MinIndex] = min(Distance);
        if (MinVal <= MinDistance) & (abs(WormSizes(MinIndex) - Tracks(ActiveTracks(i)).LastSize) < SizeChangeThreshold)
            Tracks(ActiveTracks(i)).Path = [Tracks(ActiveTracks(i)).Path; WormCoordinates(MinIndex, :)];
            Tracks(ActiveTracks(i)).LastCoordinates = WormCoordinates(MinIndex, :);
            Tracks(ActiveTracks(i)).Frames = [Tracks(ActiveTracks(i)).Frames, Frame];
            Tracks(ActiveTracks(i)).Size = [Tracks(ActiveTracks(i)).Size, WormSizes(MinIndex)];
            Tracks(ActiveTracks(i)).LastSize = WormSizes(MinIndex);
            Tracks(ActiveTracks(i)).Eccentricity = [Tracks(ActiveTracks(i)).Eccentricity, WormEccentricities(MinIndex)];
            Tracks(ActiveTracks(i)).MajorAxes = [Tracks(ActiveTracks(i)).MajorAxes, WormMajorAxes(MinIndex)];
            Tracks(ActiveTracks(i)).MinorAxes = [Tracks(ActiveTracks(i)).MinorAxes, WormMinorAxes(MinIndex)];
            Tracks(ActiveTracks(i)).Orientation = [Tracks(ActiveTracks(i)).Orientation, WormOrientation(MinIndex)];
            Tracks(ActiveTracks(i)).Box = [Tracks(ActiveTracks(i)).Box; WormBoundingBox(MinIndex,:)];
            TrackFrameNum = length(Tracks(ActiveTracks(i)).Size);
            Tracks(ActiveTracks(i)).Frame(TrackFrameNum).Image = WormImage(MinIndex).Image;

            WormCoordinates(MinIndex,:) = NaN;
        else
            Tracks(ActiveTracks(i)).Active = 0;
            if length(Tracks(ActiveTracks(i)).Frames) < MinTrackLength
                Tracks(ActiveTracks(i)) = [];
                ActiveTracks = ActiveTracks - 1;
            end
        end
    end

    % Start new tracks for coordinates not assigned to existing tracks
    NumTracks = length(Tracks);
    for i = 1:length(WormCoordinates(:,1))
        Index = NumTracks + i;
        Tracks(Index).Active = 1;
        Tracks(Index).Path = WormCoordinates(i,:);
        Tracks(Index).LastCoordinates = WormCoordinates(i,:);
        Tracks(Index).Frames = Frame;
        Tracks(Index).Size = WormSizes(i);
        Tracks(Index).LastSize = WormSizes(i);
        Tracks(Index).Eccentricity = WormEccentricities(i);
        Tracks(Index).MajorAxes = WormMajorAxes(i);
        Tracks(Index).MinorAxes = WormMinorAxes(i);
        Tracks(Index).Orientation = WormOrientation(i);
        Tracks(Index).Box = WormBoundingBox(i,:);
        Tracks(Index).Frame(1).Image = WormImage(i).Image;
    end

    % Display every PlotFrameRate'th frame
    if ~mod(Frame, PlotFrameRate)
        t1 = toc;
        PlotFrame(WTFigH, Mov, Tracks);

        FigureName = ['Movie ',num2str(MovieNum),': ',File(MovieNum).MovieFile, ...
            ' - Frame ', num2str(Frame)];
        set(WTFigH, 'Name', FigureName);

        [Level,AnimalPix] = AutoThreshold(Movdiv);
        Settings = [Settings; Frame, Level, AnimalPix];

        t2 = toc; tic;
        fps = PlotFrameRate/t1;
        fprintf('\nFrame: %5d - Time: %1.3f fps (%1.2f s) Level: (%1.3f/%3d)',Frame,fps,t2-t1,Level,AnimalPix)
         
        if (t1/PlotFrameRate) > 10     % stop if it takes too long to analyze- probably an error
            break;
        end
    end

    if TrackDye
        Movdata = Mov(:,:,1);
        bgimage = imcrop(Movdata,File(MovieNum).BgBox);
        bgavg = mean(reshape(bgimage,[],1));

        for a = 1:NumArenas
            upimage = imcrop(Movdata,File(MovieNum).Arena(a).UpBox);
            downimage = imcrop(Movdata,File(MovieNum).Arena(a).DnBox);

            uplist = sort(reshape(upimage,[],1),'descend');
            upavg = mean(uplist(1:round(length(uplist)/2))); %just average the brightest 1/2 of pixels (in case a worm is in field)

            downlist = sort(reshape(downimage,[],1),'descend');
            downavg = mean(downlist(1:round(length(downlist)/2)));

            DyeData.Arena(a).Up = [DyeData.Arena(a).Up, upavg-bgavg];
            DyeData.Arena(a).Down = [DyeData.Arena(a).Down, downavg-bgavg];
        end

        DyeData.Time = [DyeData.Time, Frame/File(MovieNum).FrameRate];
        DyeData.Background = [DyeData.Background, bgavg];
    end

    if mod(Frame,6) == 0 fprintf('-'); end
    if mod(Frame,6) == 1 fprintf('\b\\'); end
    if mod(Frame,6) == 2 fprintf('\b|'); end
    if mod(Frame,6) == 3 fprintf('\b/'); end
    if mod(Frame,6) == 4 fprintf('\b='); end
    if mod(Frame,6) == 5 fprintf('\b+'); end

end

% Get rid of invalid tracks
DeleteTracks = [];
for i = 1:length(Tracks)
    if length(Tracks(i).Frames) < MinTrackLength
        DeleteTracks = [DeleteTracks, i];
    end
end
Tracks(DeleteTracks) = [];

% Save Tracks
for a = 1:NumArenas
    DyeData.Arena(a).UpBox = File(MovieNum).Arena(a).UpBox;
    DyeData.Arena(a).DnBox = File(MovieNum).Arena(a).DnBox;
    DyeData.Arena(a).TrackBox = File(MovieNum).Arena(a).TrackBox;
end
DyeData.BgBox = File(MovieNum).BgBox;

ExpData.PixelSize = File(MovieNum).PixelSize;
ExpData.ArenaSize = File(MovieNum).ArenaSize;
ExpData.FrameRate = File(MovieNum).FrameRate;
ExpData.MovieNum = MovieNum;
ExpData.TrackTime = datestr(now);
ExpData.TrackedFrames = File(MovieNum).TrackFrames - File(MovieNum).StartFrame + 1;
ExpData.TrackStats.MinWormArea = MinWormArea;  
ExpData.TrackStats.MaxWormArea = MaxWormArea;
ExpData.TrackStats.AnimalPix = AnimalPix;  
ExpData.TrackStats.Level = Level;
ExpData.TrackStats.MinDistance = MinDistance;
ExpData.TrackStats.SizeChangeThreshold = SizeChangeThreshold;
ExpData.TrackStats.MinTrackLength = MinTrackLength;
ExpData.TrackSettings = Settings;

% Save Fragment File
disp([datestr(now),' Saving Data for Movie ',num2str(MovieNum),', fragment ',num2str(Fragment)]);
save(char(FragmentSaveNames(Fragment)), 'Tracks', 'AllData', 'DyeData', 'background', 'ExpData', 'File');
disp([datestr(now),' *** Save complete *** ']);
