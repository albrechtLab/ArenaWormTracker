%
% Combines track fragments together into one file.
%
% USAGE:
%   result = stitchtracks(inputfiles,outputfile)
%
%   inputfiles: fragment filenames
%   outputfile: output data filenames

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 29-Mar-2011 22:49:00 
%---------------------------- 

function result = stitchtracks(inputfiles,outputfile)

% IF no input arguments...
if nargin < 1

    inputfiles = [];
    [FileName, PathName] = uigetfile('*.mat', 'Select Track Files For Analysis','MultiSelect','on');
    if length(FileName) == 0 & ~exist('Tracks')
        errordlg('No file was selected for analysis');
        return;
    end
    cd(PathName);
    for fn = 1:length(FileName);
        FullName = [PathName, FileName{fn}];
        inputfiles = [inputfiles; {FullName}]; 
    end
end

if nargin < 2
    % Get output save name 
    % --------------------

    [SaveName, PathName] = uiputfile('*.mat', 'Save Track file as', [char(inputfiles(1)),'_out']);  
    outputfile = [PathName, SaveName];
end


numfiles = length(inputfiles);

result.Input = inputfiles;
result.Output = outputfile;
result.numfiles = numfiles;

if numfiles == 1 
    copyfile(char(inputfiles(1)),char(outputfile));
else
    for filenum = 1:numfiles
        disp([datestr(now),' -- Loading Track File ',int2str(filenum),'... ']);
            load(char(inputfiles(filenum)),'AllData','Tracks','background','DyeData','ExpData','File');
        disp([datestr(now),' --  ...data loaded.']);

        if filenum == 1
            StitchAllData = AllData;
            StitchTracks = Tracks;
            StitchBG = background;
            StitchDyeData = DyeData;
            StitchExpData = ExpData;
            StitchFile = File;
        else
            StitchAllData = [StitchAllData; AllData];
            StitchTracks = [StitchTracks, Tracks];
            StitchDyeData.Time = [StitchDyeData.Time, DyeData.Time];
            StitchDyeData.Background = [StitchDyeData.Background, DyeData.Background];
            for a = 1:File(ExpData.MovieNum).NumArenas
                StitchDyeData.Arena(a).Up = [StitchDyeData.Arena(a).Up, DyeData.Arena(a).Up];
                StitchDyeData.Arena(a).Down = [StitchDyeData.Arena(a).Down, DyeData.Arena(a).Down];
            end
        end
    end

    AllData = StitchAllData;
    Tracks = StitchTracks;
    background = StitchBG;
    DyeData = StitchDyeData;
    ExpData = StitchExpData;
    File = StitchFile;

    %--------------------------
    % Save stitched data file 
    %--------------------------
    disp([datestr(now),' -- Saving Stitched Track Files 1-',int2str(numfiles),'... ']);
    save(char(outputfile),'AllData','Tracks','background','DyeData','ExpData','File');
    disp([datestr(now),' --  ...file saved.']);
end


    