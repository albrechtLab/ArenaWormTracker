%
% Analyzes position, direction, and type of turns.
%
% USAGE:
%   [TurnData,SegData,SegPos,SegDir,segidx] = getTurnData(Tracks)
% 
%   Tracks: structure from ArenaTracker script
%   TurnData:
%        1. time(fr)
%        2. posX(pix)
%        3. posY(pix)
%        4. turnID(behaviorCode)
%        5. priorDirUD(1=up(-y), 2=horiz, 3=down)
%        6. afterDirUD(1=up(-y), 2=horiz, 3=down)
%   SegData:
%        1. Track #
%        2. Segment # (of Track)
%        3. Behavior Code 
%        4. Behavior Run (1=Run, 0=Turn)
%        5. Segment Length (frames)
%   SegPos:
%      1-3. SegmentStart (frame, X, Y of starting point)
%      4-6. SegmentEnd (frame, X, Y of ending point)
%      7-8. SegmentCentroid (X, Y of centroid)
%     9-12. SegmentMinMax (minX, minY, maxX, maxY)
%   SegDir
%        1. EndAngleUD(1=up(-y), 2=horiz, 3=down)
%      2-3. SegmentDir (mean, std of all angles in segment: +x axis; ccw)
%        4. SegmentBoundHead (initial, final angle)
%     5-17. Histogram of angles (# frames in 0,30,60,90,... 330)
%   segidx(turnidx) = segment # of forward run prior to turn

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 04-Oct-2009 23:44:27 
%---------------------------- 

function [TurnData,SegData,SegPos,SegDir,segidx] = getTurnData(Tracks)

%----------------------------------------------%
%               SEGMENT ANALYSIS               %
%----------------------------------------------%
tic;

SegData = []; SegPos = []; SegDir = []; 
anglegroup = [2 1 1 2 3 3];  % grouping for 0-60-120-180-240-300 degrees
endstallframes = 10; % max number of frames to test for stalling

NumTracks = length(Tracks);
for tr = 1:NumTracks
    StartFrame = Tracks(tr).Frames(1);
    Segments = length(Tracks(tr).Segment);
    SegmentTrack = tr*ones(Segments,1);
    SegmentCode = [Tracks(tr).Segment.Code]';
    SegmentRun = [Tracks(tr).Segment.Run]';
    SegmentLength = [Tracks(tr).Segment.Length]';
    SegmentSeg = (1:Segments)';
    
    %analyze stalls
    nostall = Tracks(tr).Frames; stall = Tracks(tr).Stall; nostall(stall) = NaN; nostall = sort(nostall); nostall = nostall(1:length(nostall)-length(stall));
    [a,b] = findsegment(stall);
    frSeg = []; for s=1:Segments; frSeg = [frSeg,ones(1,SegmentLength(s))*s]; end
    trackframe = [];
    trackframe(nostall - StartFrame + 1) = frSeg; trackframe(stall) = 0;
    tracklength = length(trackframe);
    c = [];
    for st = 1:size(a,2)
        st1 = a(st)-b(st)+1; st2 = a(st);
        segbef = NaN; segaft = NaN;
        if st1 > 1              segbef = trackframe(st1 - 1); end
        if st2 < tracklength    segaft = trackframe(st2 + 1); end
        c = [c, nanmean([segbef, segaft])];
    end
    StallData = [a',b',c'];
    segfr = zeros(Segments,1); for s=1:Segments; segfr(s)=max(find(trackframe == s)); end
    segfr = [[1; segfr(1:Segments-1)+1], segfr];
    SegmentEndStall = zeros(Segments,1);
    for s = 1:Segments
        segstall = (trackframe(max(segfr(s,1),segfr(s,2)-endstallframes+1):segfr(s,2)) == 0);
        SegmentEndStall(s) = sum(segstall .* 2.^(length(segstall)-1:-1:0));
    end
    
    SegmentStall = [];
    SegmentStart = []; SegmentEnd = []; SegmentCentroid = []; SegmentDir = []; 
    SegmentAngHist = []; SegmentMinMax = []; SegmentBoundHead = [];
    for seg = 1:Segments
        if size(StallData,1) > 0
            SegmentStall = [SegmentStall; sum(StallData(find(StallData(:,3) == seg-0.5),2)), ...
                sum(StallData(find(StallData(:,3) == seg),2)), sum(StallData(find(StallData(:,3) == seg+0.5),2))];
        else
            SegmentStall = zeros(Segments,3);
        end
        SegmentStart = [SegmentStart; nostall(Tracks(tr).Segment(seg).Frames(1)), Tracks(tr).Segment(seg).Path(1,:)];
        endfr = Tracks(tr).Segment(seg).Length;
        SegmentEnd = [SegmentEnd; nostall(Tracks(tr).Segment(seg).Frames(endfr)), Tracks(tr).Segment(seg).Path(endfr,:)];
        SegmentCentroid = [SegmentCentroid; Tracks(tr).Segment(seg).Centroid];
        SegmentMinMax = [SegmentMinMax; min(Tracks(tr).Segment(seg).Path,[],1), max(Tracks(tr).Segment(seg).Path,[],1)];
        
        ForwardDir1 = rem(Tracks(tr).Segment(seg).PathAngle' + ~(Tracks(tr).Segment(seg).Forward')*180,360);
        ForwardDir2 = rem(Tracks(tr).Segment(seg).PathAngle' + ~(Tracks(tr).Segment(seg).Forward')*180 + 180,360)-180;
        std1 = std(ForwardDir1); std2 = std(ForwardDir2);
        if std1 < std2  % use option with lower stdev
            SegmentDir = [SegmentDir; mean(ForwardDir1), std(ForwardDir1)];
        else
            SegmentDir = [SegmentDir; mean(ForwardDir2), std(ForwardDir2)];
        end
        
        anghist = hist(Tracks(tr).Segment(seg).PathAngle,0:30:360); anghist = [anghist(1)+anghist(13), anghist(2:12)];
        SegmentAngHist = [SegmentAngHist; anghist];
        
        boundangle = Tracks(tr).Segment(seg).BoundAngle;
        boundforward = Tracks(tr).Segment(seg).BoundForward;
        boundheading = rem(boundangle + 180*(~boundforward),360);
        SegmentBoundHead = [SegmentBoundHead; boundheading];
    end
    
    angUD = anglegroup((floor(SegmentBoundHead(:,2)/60)+1))' .* SegmentRun;
    
    SegData = [SegData; [SegmentTrack, SegmentSeg, SegmentCode, SegmentRun, ...
                         SegmentLength, SegmentStall, SegmentEndStall]];
    SegPos = [SegPos; [SegmentStart, SegmentEnd, SegmentCentroid, SegmentMinMax]];
    SegDir = [SegDir; [angUD, SegmentDir, SegmentBoundHead, SegmentAngHist]];
end

runseg = SegData(:,4) & SegDir(:,4)==SegDir(:,5);

t = toc; disp(['Segment analysis: ',num2str(t),' s']);
tic;

nseg = size(SegPos,1);
f0 = 1:nseg-2;  % forward segment BEFORE turn 
f1 = 2:nseg-1;  % turn
f2 = 3:nseg;    % forward segment AFTER turn 
    
segidx = find(runseg(f0)==1 & runseg(f2)==1 & SegData(f0,1) == SegData(f2,1) & ...  % pre and post segments are forward runs of the same track
        SegData(f1,3) ~= 2);   

TurnData = [];
for tnum = 1:size(segidx,1)
    TurnData = [TurnData; SegPos(segidx(tnum),4:6), SegData(segidx(tnum)+1,3), ...  % 1. time(fr), 2. posX, 3. posY, 4. turnID
        SegDir(segidx(tnum),1), SegDir(segidx(tnum)+2,1), ...                       % 5. priorDirUD, 6. afterDirUD
        sum(SegData(segidx(tnum),6:7)), sum(SegData(segidx(tnum)+1,6:7))];          % 7. priorStallFr, 8. turnStallFr
end

t = toc; disp(['Turn analysis: ',num2str(t),' s']);


