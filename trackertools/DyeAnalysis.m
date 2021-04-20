%
% Analyze dye concentration to determine flow properties in a pulse assay.
%
% USAGE:
%   output = DyeAnalysis(DyeData,ontime)
%
%   DyeData: dye concentration data from ArenaTracker script
%   ontime: n x 2 vector contining dye on-off frame numbers
%           e.g. [10 20; 30 40] for dye on at frame 10; off at 20; on at 30; etc.

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 30-Mar-2011 14:29:13 
%---------------------------- 

function output = DyeAnalysis(DyeData,ontime)

if ~exist('ontime') || isempty(ontime) ontime = NaN; end

numarenas = length(DyeData.Arena);

for a = 1:numarenas
    xup = DyeData.Arena(a).UpBox(1)+DyeData.Arena(a).UpBox(3)/2;
    xdown = DyeData.Arena(a).DnBox(1)+DyeData.Arena(a).DnBox(3)/2;
    up = DyeData.Arena(a).Up; dn = DyeData.Arena(a).Down;

    DyePresent = (max(smooth(up)) - min(smooth(up))) > 8; % changed from 15 20130909 DRA

    if DyePresent & ~isnan(ontime) % added 20140625 check that ontime exists for pulse assay
        dyelen = length(DyeData.Time);
        updye = smooth(scale(up,0,1))';
        downdye = smooth(scale(dn,0,1))';

        upswitch = find(sign(updye(1:dyelen-1)-0.5) ~= sign(updye(2:dyelen)-0.5));
        downswitch = find(sign(downdye(1:dyelen-1)-0.5) ~= sign(downdye(2:dyelen)-0.5));

        if upswitch(1) > downswitch(1) downswitch = circshift(downswitch,[0 -1]); end

        switchlen = [length(upswitch), length(downswitch)];
        if switchlen(1) ~= switchlen(2)
            upswitch = upswitch(1:min(switchlen));
            downswitch = downswitch(1:min(switchlen));
        end
        timedelay = [upswitch; downswitch; downswitch-upswitch];
        avgdelay = mean(timedelay(3,:));  %upstream to downstream delay (fr)
        stddelay = std(timedelay(3,:));

        if ~isnan(ontime)
            startdelay = upswitch - sort(ontime(1:length(upswitch)));
            avgstartdelay = mean(startdelay); % delay from switch to upstream point (fr)
            stdstartdelay = std(startdelay);

            if (stdstartdelay > avgstartdelay)  % this happens when StreamPix loses frames
                avgstartdelay = mean(startdelay(1:(length(startdelay)/2))); % only use first 1/2 of time
                stdstartdelay = std(startdelay(1:(length(startdelay)/2)));
                frameslost = '***';
            end
            output(a).StartDelayFr = [avgstartdelay, stdstartdelay];
        else
            output(a).StartDelayFr = NaN;
        end
       
        flowvelocity = abs(xdown - xup) ./ timedelay(3,:); % pix/fr linear velocity
        avgflow = mean(flowvelocity);
        stdflow = std(flowvelocity);

        output(a).TimeDelayFr = timedelay;
        output(a).UpDownDelayFr = [avgdelay, stddelay];
        output(a).VelocityPixPerFr = [avgflow, stdflow];
    else
        output(a).StartDelayFr = NaN;
        output(a).TimeDelayFr = NaN;
        output(a).UpDownDelayFr = NaN;
        output(a).VelocityPixPerFr = NaN;
    end

    output(a).XUpDownPix = [xup, xdown];

%         disp(['*** Dye Front Tracked ***', datestr(now)]); 
end    
