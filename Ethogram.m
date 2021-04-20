%
% Displays ethogram of behavior and summarizes behavioral state probability
% and speed over time.  Adjusts timing according to flow properties, either
% via user input or automatically from dye experiments. Saves data to a
% file (*_ethogram.mat) and prints PDF summary pages.    
%
% USAGE:
%   Ethogram(FileName)
%
%   FileName: single filename of segmented data or cell array of multiple
%               filenames. Select with user input if none given.

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 31-Mar-2011 12:09:13 
%---------------------------- 

function Ethogram(FileName)
    
if nargin < 1 FileName = {}; end
if length(FileName) == 0
    % Get track data for analysis
    % --------------------------
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

        SaveName = strrepl(FullName,'_seg.mat','_ethogram.mat');
        PulseProgramName = strrepl(FullName,'_seg.mat','.txt');

        %------------
        % Load Data
        %------------
        clear('Tracks','DyeData','ExpData','numcycles');

        disp([datestr(now),': Loading File #',num2str(fnum),': ',FullName]);

        load(FullName,'Tracks','DyeData','ExpData');

        try
            stimulusdata = parseExptData(PulseProgramName);
        catch
            stimulusdata = [];
            disp(['** No experiment data file found: ',PulseProgramName]);
            disp('Results will not be time-shifted');
        end

        %--------------------
        % Analyze Dye Front
        %--------------------

        if ~isempty(stimulusdata)
            ontime = struct2mat(1,stimulusdata.Stimulus,[],{'Time'});
            if isfield(stimulusdata,'NumCycles') 
                numcycles = str2num(char(stimulusdata.NumCycles));
                if isfield(stimulusdata,'CycleLen') 
                    cyclelen = str2num(char(stimulusdata.CycleLen)); 
                else
                    cyclelen = max(ontime);
                end
                cycleontime = ontime;
                for i=1:numcycles-1;
                    ontime = [ontime; cycleontime+i*cyclelen];
                end
            end
        else
            ontime = NaN;
        end

        flowinfo = DyeAnalysis(DyeData,ontime*60*ExpData.FrameRate);
        ExpData.Flow = flowinfo(ExpData.Arena);
        ExpData.Flow.VelocityMMperSec = ExpData.Flow.VelocityPixPerFr * ExpData.FrameRate / ExpData.PixelSize;
        ExpData.Flow.StartDelaySec = ExpData.Flow.StartDelayFr / ExpData.FrameRate;

        % use velocity and start delay data from experiment info file
        if isfield(stimulusdata,'FlowVel') 
            ExpData.Flow.VelocityMMperSec = str2num(char(stimulusdata.FlowVel));
            ExpData.Flow.VelocityPixPerFr = ExpData.Flow.VelocityMMperSec * ExpData.PixelSize / ExpData.FrameRate; 
        end
        if isfield(stimulusdata,'Delay')
            ExpData.Flow.StartDelaySec = str2num(char(stimulusdata.Delay));
            ExpData.Flow.StartDelayFr = ExpData.Flow.StartDelaySec * ExpData.FrameRate;
        end
            
        %--------------------------------------
        % Generate ethogram and data matrices
        %--------------------------------------

        ExpData.Flow.DelayFr = ExpData.Flow.StartDelayFr(1) - ExpData.Flow.XUpDownPix(1) / ExpData.Flow.VelocityPixPerFr(1);

        ethfig = findobj(get(0,'Children'),'Tag','Ethogram');
        if isempty(ethfig) ethfig = figure; set(ethfig,'Tag','Ethogram'); end

        Data = Tracks2Matrix(Tracks,ExpData.Flow,~isnan(ExpData.Flow.DelayFr));

        %------------------------------
        % Save Data and summary plots
        %------------------------------

        save(SaveName,'Data','ExpData');
        
        % page 1: Raw matrices
        orient(ethfig,'tall');
        saveas(ethfig,strrepl(SaveName,'.mat','_raw.pdf'),'pdf');

        % page 2: Population average state probability and speed
        cmap = [1 1 1;.7 .7 .7; .7 .7 .7; 0 0 0; .3 .3 .3; .6 0 0; 1 .2 .2; 1 1 1; .9 .9 .9];
        clf;
        t = (1:size(Data.behmat,2)) / 60 / ExpData.FrameRate; % min
        subplot(5,1,1); image(t,1:size(Data.behmat,1),ind2rgb(Data.behmat,cmap)); ylabel('Animal #');
        if isnan(ExpData.Flow.VelocityMMperSec) 
            FlowLabel = '(real time)';
        else
            FlowLabel = sprintf('(flow: %1.3f mm/s, delay: %1.3f s)',ExpData.Flow.VelocityMMperSec(1),ExpData.Flow.StartDelaySec(1));
        end
        title(sprintf('%s %s',SaveName,FlowLabel),'interpreter','none'); 
        subplot(5,1,2); PlotStates(Data.behprob,[],t,0,0,0); ylabel('State probability');
        subplot(5,1,[3:4]); PlotStates(Data.behprob,[],t,.05,0,0); hilite(ontime,[],[1 1 .5]); ylabel('State probability'); xlim([0 t(end)]);
        subplot(5,1,5); plot(t,Data.speed.fwdpause); ylim([0 0.4]); hilite(ontime,[],[1 1 .5]); ylabel('Speed (mm/s)'); xlim([0 t(end)]);
        xlabel('Time (min)');
        orient(ethfig,'tall');
        saveas(ethfig,strrepl(SaveName,'.mat','_avg.pdf'),'pdf');

        % page 3: Population state probability and speed, averaged over
        % repeated cycles
        if exist('numcycles') && numcycles > 1
            clf;
            cycbehmat = split(Data.behmat,4,2,1);
            t = (1:size(cycbehmat,2)) / 60 / ExpData.FrameRate; % min
            subplot(5,1,1); image(t,1:size(cycbehmat,1),ind2rgb(cycbehmat,cmap)); ylabel('Animal #');
            title(sprintf('%s %s',SaveName,FlowLabel),'interpreter','none'); 
            cycbehhist = hist(cycbehmat,1:8);
            cycbehnum = sum(cycbehhist(1:6,:));
            cycbehprob = cycbehhist(1:6,:) ./ repmat(cycbehnum,6,1); 
            subplot(5,1,2); PlotStates(cycbehprob,[],t,0,0,0); ylabel('State probability');
            subplot(5,1,[3:4]); PlotStates(cycbehprob,[],t,.05,0,0); hilite(ontime,[],[1 1 .5]); ylabel('State probability');  xlim([0 t(end)]);
            subplot(5,1,5); plot(t,nanmean(reshape(Data.speed.fwdpause',[],numcycles)')); ylim([0 0.4]); hilite(ontime,[],[1 1 .5]); ylabel('Speed (mm/s)'); xlim([0 t(end)]); 
            xlabel('Time (min)');
            orient(ethfig,'tall');
            saveas(ethfig,strrepl(SaveName,'.mat','_cycavg.pdf'),'pdf');
        end
end
