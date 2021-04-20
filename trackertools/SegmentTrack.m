%
% Segment the specified track according to behavior.
%
% USAGE:
%   [Segment, SegAnalysis] = SegmentTrack(CurrentTrack,Settings)
%
%   CurrentTrack: current element of Tracks structure
%   Settings: segmentation settings structure
%   Segment: structure containing data for each track segment
%   SegAnalysis: structure of post-segmentation analysis

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 10-Feb-2011 16:34:55 
%---------------------------- 

function [Segment, SegAnalysis] = SegmentTrack(CurrentTrack,Settings)

if nargin < 2
    SegmentationSettings;
end

se1 = strel('line',5,0); se2 = strel('line',3,0); % structure elements to open/close gaps of 1-2 or <=4

%tic;
xy = CurrentTrack.Path;
frames = CurrentTrack.Frames;
trlength = length(frames);
bodyangle = mod([CurrentTrack.Orientation]+30,60)-30;      % make the range from -30 to +30 degrees
bodyangle360 = mod([CurrentTrack.Orientation]+30,360)-30;  % make the range from -30 to +330 degrees
ecc = sqrt(1 - CurrentTrack.Eccentricity.^2);              % eccentricity = minor axis / major axis
mindistance = Settings.StallDistance * Settings.PixelSize / Settings.FrameRate;  % max pixel movement per frame for stall

[smang, smangvel, smangacc,   smdist] = angle(xy);   % based on smoothed centroid data as mean between adjacent
[  ang,   angvel,   angacc, origdist] = angle2(xy);  % based on raw centroid data
niceangle = ang + 360*(ang < 0);        % make the range from   0 to +360 degrees
niceangle = mod(niceangle+30,360)-30;   % make the range from -30 to +330 degrees
pathangle = mod(ang+30,60)-30;          % make the range from -30 to +30 degrees
smpathangle = mod(smang+30,60)-30;

%----------------------
% Find stalled frames
%----------------------
stallfr = find(origdist > 0 & origdist <= mindistance);
if (length(stallfr) > 0 & Settings.FixJitter)
    
    %get rid of stalled frames and recalculate variables
    [a,nostallfr] = veccomp(1:trlength,stallfr);
    xy = xy(nostallfr,:);
    frames = frames(nostallfr);
    trlength = length(frames);
    bodyangle = bodyangle(nostallfr);
    bodyangle360 = bodyangle360(nostallfr);
    ecc = ecc(nostallfr);

    [smang, smangvel, smangacc, smdist] = angle(xy);
    [  ang,   angvel,   angacc,   dist] = angle2(xy);
    niceangle = ang + 360*(ang < 0);        % make the range from   0 to +360 degrees
    niceangle = mod(niceangle+30,360)-30;   % make the range from -30 to +330 degrees
    pathangle = mod(ang+30,60)-30;          % make the range from -30 to +30 degrees
    smpathangle = mod(smang+30,60)-30;
   
else
    nostallfr = 1:trlength;
    dist = origdist;
end
    

windowSize = Settings.SmoothWinSize;
fpathangle = filter(ones(1,windowSize)/windowSize,1,smpathangle');
fbodyangle = filter(ones(1,windowSize)/windowSize,1,bodyangle);

% Analyze animal size
relativesize = CurrentTrack.Size / median(CurrentTrack.Size);
collision = relativesize > Settings.CollisionRelSize;

relativelength = CurrentTrack.MajorAxes / max(CurrentTrack.MajorAxes(~collision));

trackbreakfr = [];
if (trlength >= 10) & 0
    
    scrunching = imclose(relativelength <= 0.6,strel('line',3,0));
    [a,b] = findsegment(find(scrunching));
    trackbreakfr = a-b+8;
    trackbreakfr = trackbreakfr(find(b >= 2)); % in real frames
    
    if length(trackbreakfr) > 0
        neartostall = min(abs(addortho(stallfr,-trackbreakfr)),[],1); % find distance to nearest stalled frame
        validbreakpoint = find(neartostall > 10);
        if length(validbreakpoint) > 0
            trackbreakfr = trackbreakfr(validbreakpoint);

            neartonostall = abs(addortho(nostallfr',-trackbreakfr)); % find nearest non-stalled frame
            [a,b] = sort(neartonostall);
            trackbreakfr = b(1,:);
        else
            trackbreakfr = [];
        end
    end
end


%--------------------------------------------------
% Identify segments based on straight forward runsstallfr
%--------------------------------------------------

strindex = find(abs(fpathangle-fbodyangle) < Settings.MaxBodyPathAngleDev & ...
                abs(fpathangle) < Settings.MaxPathAngleDev & ...
                ecc < Settings.MaxFwdRunEcc);

[a,b,strindexlen] = findsegment(strindex);
longindex = strindex(find(strindexlen >= Settings.MinFwdRunFr));   % need consecutive straight frames for run

longindex = sort(unique([longindex, trackbreakfr])); % divide up tracks by straight segments OR breakpoints

if length(longindex) > 0
    [a,b,strindexlen] = findsegment(longindex);
    c = a-b+1; shift = c(1)>1;
    lena = length(a);

    segments = ones(1,trlength)*2*lena + shift;
    if shift segments(1:c(1)-1) = 1; end
    for i = 1:lena-1
        segments(c(i):a(i)) = 2*i - 1 + shift;
        segments(a(i)+1:c(i+1)-1) = 2*i + shift;
    end
    segments(c(lena):a(lena)) = 2*lena - 1 + shift;
    numsegments = max(segments);
    cmap = jet(numsegments);

    for seg = 1:numsegments
        seglist = find(segments == seg);
        straightseg(seg) = length(find(longindex == seglist(1)));
    end
else
    segments = ones(1,trlength);
    numsegments = 1;
    straightseg = 0;
end

realfrsegments = zeros(size(CurrentTrack.Frames));
realfrsegments(nostallfr) = segments;
realfrsegments(1) = 1;
for i=stallfr
    if realfrsegments(i)==0 realfrsegments(i)=realfrsegments(i-1); end
end


%--------------------- 
% collect raw data
%---------------------

for seg = 1:numsegments
    
    Segment(seg).Frames = find(segments == seg);
    Segment(seg).RealFrames = nostallfr(Segment(seg).Frames);
    [a, Segment(seg).RealStallFrames] = veccomp(find(realfrsegments == seg),Segment(seg).RealFrames);
    Segment(seg).Path = xy(Segment(seg).Frames,:);
    Segment(seg).Centroid = mean(Segment(seg).Path,1);
    Segment(seg).Length = length(Segment(seg).Frames);
    
    Segment(seg).Run = straightseg(seg) + 2*(Segment(seg).Length == 1);    % classify forward straight runs (2 = trackbreak)
    
    Segment(seg).Distance = dist(Segment(seg).Frames)';
    Segment(seg).Orientation = bodyangle360(Segment(seg).Frames);
    Segment(seg).PathAngle = niceangle(Segment(seg).Frames)';
    Segment(seg).PathAngVel = angvel(Segment(seg).Frames)';

    Segment(seg).PathAngAcc = angacc(Segment(seg).Frames)';
    
    if Segment(seg).Run
        Segment(seg).Direction = round(Segment(seg).PathAngle / 60) * 60;   % forward runs have 60' directions
    else
        Segment(seg).Direction = round(Segment(seg).PathAngle / 30) * 30;   % turns may not...
        Segment(seg).Direction = Segment(seg).Direction + 360 * (Segment(seg).Direction < 0);
    end
    Segment(seg).MeanAngle = mean(Segment(seg).PathAngle');
    Segment(seg).StdAngle = std(Segment(seg).PathAngle');
    
    Segment(seg).BoundAngle = [Segment(seg).Direction(1), Segment(seg).Direction(Segment(seg).Length)];
    if seg > 1
        if ~Segment(seg).Run Segment(seg).BoundAngle(1) = Segment(seg-1).Direction(Segment(seg-1).Length); end
        if ~Segment(seg-1).Run Segment(seg-1).BoundAngle(2) = Segment(seg).Direction(1); end
    end
    
    if Segment(seg).Length > 1
        Segment(seg).SmoothAngChange = (0:1/(Segment(seg).Length - 1):1) * ...
                    (Segment(seg).PathAngle(Segment(seg).Length) - Segment(seg).PathAngle(1)) + ...
                     Segment(seg).PathAngle(1);
    else
        Segment(seg).SmoothAngChange = Segment(seg).PathAngle(1);
    end
    Segment(seg).PathAngle = Segment(seg).PathAngle + 360 * (Segment(seg).PathAngle < 0);   
    Segment(seg).Collision = any(collision(Segment(seg).RealFrames));
    Segment(seg).MajorAxes = CurrentTrack.MajorAxes(Segment(seg).RealFrames) / ...
        (prctile(CurrentTrack.MajorAxes .* ~imdilate(collision,strel('line',15,0)),99)+eps);
end

%------------------------------------------
% analyze straight forward 'Run' segments
%------------------------------------------
forwardturn = []; 

runseg = find([Segment.Run]);
for seg = runseg
    
    Segment(seg).TrackDirection = mod(Segment(seg).Direction,180);
    Segment(seg).Forward = (Segment(seg).TrackDirection == Segment(seg).Direction);
    
    if abs(Segment(seg).PathAngVel(end)) >= (Settings.MinRevAngVel / Settings.FrameRate)
        Segment(seg).Forward(end) = ~Segment(seg).Forward(end);
    end
      
    Segment(seg).DirectionChange = mod([0, Segment(seg).Direction(2:Segment(seg).Length) - Segment(seg).Direction(1:Segment(seg).Length-1)]+180,360)-180;

    turnfr = find(abs(Segment(seg).DirectionChange) == 60);
    swapprior = find(Segment(seg).Forward(turnfr) ~= Segment(seg).Forward(turnfr-1) & ...
                     abs(sum(Segment(seg).DirectionChange(turnfr:min(Segment(seg).Length,turnfr+1)))) < 180);
    for fturn = 1:length(swapprior)
        Segment(seg).Forward(1:turnfr(swapprior)-1) = ~Segment(seg).Forward(1:turnfr(swapprior)-1);
    end
    
    if (sum(Segment(seg).Forward) < Segment(seg).Length * 0.55 & ...     % if more reverse than forward, flip
        (sum(Segment(seg).Forward) < Segment(seg).Length * 0.45 | Segment(seg).Forward(1) == 0))  % if equal, assume that starts out forward
        Segment(seg).Forward = ~Segment(seg).Forward;
        Segment(seg).Flip = 1;
    end
    
    % find any 60 degree turns
    Segment(seg).Turn60 = find(abs(mod((Segment(seg).Direction - Segment(seg).BoundAngle(1))+180,360)-180) == 60);
    if (length(Segment(seg).Turn60) > 2 || ...    % require >2 frames at the new angle if mid-segment
            (length(Segment(seg).Turn60) > 0 && ...  % or just 1 frame if at the beginning or end of a segment
            (min(Segment(seg).Turn60) == 1 || max(Segment(seg).Turn60) == Segment(seg).Length)))
        turnframe = Segment(seg).Frames(Segment(seg).Turn60(1));
        forwardturn = [forwardturn; turnframe, xy(turnframe,:), Segment(seg).Direction(Segment(seg).Turn60(1)) - Segment(seg).BoundAngle(1)];
    end
    
    if sum(~Segment(seg).Forward) > 0   % if there are reversals, find them    
        if (Segment(seg).Length < 7 & sum(~Segment(seg).Forward) > sum(Segment(seg).Forward)) % if more revs than forwards, it's probably flipped
            Segment(seg).Forward = ~Segment(seg).Forward;
        end
        [rfr,rlen,rindex] = findsegment(find(Segment(seg).Forward == 0));
        Segment(seg).Reversal = [rfr', rlen'];
        Segment(seg).ID = '2: Forward Run with Reversal(s)';
    else
        Segment(seg).ID = '1: Forward Run';
    end
    Segment(seg).BoundForward = [Segment(seg).Forward(1), Segment(seg).Forward(Segment(seg).Length)];
    
    Segment(seg).Code = str2num(Segment(seg).ID(1));
end

%--------------------------
% analyze turning segments
%--------------------------
turnseg = find(~[Segment.Run]);
for seg = turnseg

    Segment(seg).Scrunch = Segment(seg).MajorAxes < 0.65; % was 0.6
    omega = (ecc(Segment(seg).Frames) > Settings.MinOmegaEcc);
    Segment(seg).Omega = imclose(omega,se2);  % filter: fill gaps of 1-2
        
    % Find Forward motion frames
    Segment(seg).Forward = ones(1,Segment(seg).Length);
    if seg > 1 
        if Segment(seg-1).BoundForward(2) == 0 
            Segment(seg).Forward = ~Segment(seg).Forward;   % start turn in rev if previous track ended in rev
            Segment(seg).BoundAngle(1) = mod(Segment(seg).BoundAngle(1) - 180, 360);
        end
    end  
    
    % deal with short reversals that may not hit the full threshold
    revfr = find(abs(Segment(seg).PathAngVel) >= (Settings.MinRevAngVel / Settings.FrameRate)); %105deg
    revfr2 = find(abs(Segment(seg).PathAngVel) >= 90); % possible reversal? if angvel = 90deg/fr, check for adjacent, opposite pairs
    [a,b] = findsegment(revfr2); 
    yesrev = (b == 1 & abs(sum(Segment(seg).PathAngVel([a; a-b+1]))) > 180);
    revfr = unique([revfr, a(yesrev)]);
    yesrev1fr = (a >= 2 & (b == 2 & (sign(Segment(seg).PathAngVel(a)) ~= sign(Segment(seg).PathAngVel(a-b+1)) | ...
                                 abs(sum(Segment(seg).PathAngVel([a; a-b+1]))) > 300))); % single frame reversal
    revfr2 = []; for r1 = find(yesrev1fr) revfr2 = [revfr2, a(r1) + [-1 0]]; end  % reconstruct revfr2 var with valid revs
    revfr = unique([revfr, revfr2]);
    
    % deal with scrunch states- animals that reverse direction between 2 posts 
    [a,b] = findsegment(find(Segment(seg).Scrunch & ~Segment(seg).Omega));
    scrunchlength = diff(Segment(seg).RealFrames(min([a-b+1; a+1],length(Segment(seg).Frames))))+1;
    scrunchlabel = bwlabel(Segment(seg).Scrunch & ~Segment(seg).Omega);
    scrunchflippedalready = []; for sc = 1:length(a); scrunchflippedalready(sc) = length(find(scrunchlabel(revfr) == sc)); end
    scrunchflipfr = a(find(scrunchlength > 12 & ~scrunchflippedalready)); 
    Segment(seg).ScrunchFlip = scrunchflipfr;
    revfr = unique([revfr, scrunchflipfr]);  % reverse if scrunchlength > 12 fr... probably a u-turn?
    
    [a,b]=findsegment(find(Segment(seg).Omega)); omegalabel = bwlabel(Segment(seg).Omega);
    omegaflips = []; for oc = 1:length(a); omegaflips(oc) = length(find(omegalabel(revfr) == oc)); end
    toomanyflips = find(b >= 2 & omegaflips == 2);
    for od = 1:length(toomanyflips); flipremove = find(omegalabel(revfr) == toomanyflips(od)); revfr = revfr(setdiff(1:length(revfr),flipremove(1))); end
    
    firstfr = Segment(seg).Forward(1);
    for r = 1:length(revfr)
        Segment(seg).Forward(revfr(r):Segment(seg).Length) = mod(r+firstfr,2);
    end
    Segment(seg).BoundForward = [Segment(seg).Forward(1), Segment(seg).Forward(Segment(seg).Length)];
    
    %********************** flip forward/reverse if necessary
    
    test1 = 0; test2 = 0; nextseglength = 0; nextsegconfidence = 0;
    if seg > 1 test1 = (Segment(seg).BoundForward(1) == ~Segment(seg-1).BoundForward(2) & ~any(revfr == 1)); end
    if seg < numsegments 
        test2 = (Segment(seg).BoundForward(2) == ~Segment(seg+1).BoundForward(1)); 
        nextseglength = Segment(seg+1).Length;
        nextsegconfidence = mean(Segment(seg+1).Forward(1:min(nextseglength,4)) == Segment(seg+1).BoundForward(1));
    end

    % if both start and end directions are wrong, flip forward frames
    % OR if first or last segment, just one wrong is enough
    [a1,rblock] = findsegment(find(~Segment(seg).Forward)); if length(rblock) == 0 rblock = 0; end % determine max run of reverse and forward
    [a1,fblock] = findsegment(find(Segment(seg).Forward)); if length(fblock) == 0 fblock = 0; end
    % flip forward if it matches adjacent segments better and it won't
    % produce a long reversal block OR flip if reversal block already too big
    maxrevfr = Settings.MaxRevTime * Settings.FrameRate;
    if ((((test1 & test2) | (test2 & seg == 1) | (test1 & seg == numsegments) | ...
            (test2 & nextseglength > 4 & nextsegconfidence > 0.6)) ...
            & max(fblock) < maxrevfr) | max(rblock) > maxrevfr)
        Segment(seg).Forward = ~Segment(seg).Forward;
        Segment(seg).BoundForward = [Segment(seg).Forward(1), Segment(seg).Forward(Segment(seg).Length)];
        Segment(seg).Flip = 1; % flip reversals to improve continuity-- check validity of this
    end
    %**********************
        
    if Segment(seg).Length <= (Settings.MaxShortTurnLen * Settings.FrameRate)     % short turning segments
        
        Segment(seg).ForwardTurn = sum(abs(Segment(seg).PathAngVel(2:Segment(seg).Length-1)) > Settings.MaxSmoothAngVel) == 0;  % no sharp direction changes
        if Segment(seg).ForwardTurn
            dirchange = mod(Segment(seg).BoundAngle(2) - Segment(seg).BoundAngle(1) + 180,360) - 180;
            if dirchange == 0
                Segment(seg).ID = '3: Forward Shift';
            elseif abs(dirchange) == 180
                Segment(seg).ID = '2a: Forward Run with Reversals (or U-turn?)';
            else
                Segment(seg).ID = ['4: Forward Turn ',num2str(dirchange),'deg ccw'];
                forwardturn = [forwardturn; mean(Segment(seg).Frames), Segment(seg).Centroid, dirchange];    
            end
        else
            Segment(seg).ID = '5: Short Reversal??';
        end
    else    % long turning segments
           
        % determine pirouette/omega turns
        
        Segment(seg).Erratic = (Segment(seg).PathAngAcc > (Settings.MinErraticAngAcc / Settings.FrameRate^2));  % find erratic movements (was 20 before)
        dirchg = [0, Segment(seg).Forward(3:Segment(seg).Length) ~= Segment(seg).Forward(1:Segment(seg).Length-2), 0];
        Segment(seg).Erratic = (Segment(seg).Erratic & ~dirchg); % remove erratic movements due to direction changes
        
        Segment(seg).Erratic = imopen(imclose(Segment(seg).Erratic,se1),se2);   % filter: fill in frames within 5; remove single or double
        
        curl = (ecc(Segment(seg).Frames) > Settings.MinPirEcc);
        Segment(seg).Curl = imopen(imclose(curl,se2),se2);  % filter: fill gaps of 1-2; remove single or double
        
        %------------
        omega2post = find(ecc(Segment(seg).Frames) >= Settings.Omega2PostMinEcc & ...
                          ecc(Segment(seg).Frames) <= Settings.MinOmegaEcc); % test for omegas around 2 posts
        if length(omega2post) > 0
            omega2postfr = min(addortho(omega2post',1:Settings.Omega2PostTestTime*Settings.FrameRate)',Segment(seg).Length);
            omega2postAngVel = Segment(seg).PathAngVel(omega2postfr);
            omega2postAngVel(find(omega2postfr > Segment(seg).Length)) = NaN;
            omega2postAngmean = nanmean(omega2postAngVel);
            Segment(seg).Omega(omega2post(find(abs(omega2postAngmean) > Settings.Omega2PostMinAngVel/Settings.FrameRate))) = 1;  % look for consistent curve
        end
        %------------
            
        Segment(seg).Pirouette = (sum(Segment(seg).Omega) >= 1) | (sum(Segment(seg).Scrunch) >= 1);          
        if Segment(seg).Pirouette
            
            % correct forward/reverse directions
            % (1) if very straight, usually reverse
            straightfr = find(ecc(Segment(seg).Frames) < 0.12);
            if length(straightfr)==0 test1 = NaN;
            else test1 = mean(Segment(seg).Forward(straightfr)) < 0.5; end
            
            % (2) after omegas, usually forward. check 3-5 frames
            [oend,olen] = findsegment(find(Segment(seg).Omega));
            if length(oend) > 0
                postomegafr = addortho(oend,(3:5)');  
            else postomegafr = []; end
            postomegafr = postomegafr(find(postomegafr <= length(Segment(seg).Frames)));
            if length(postomegafr)==0 test2 = NaN; 
            else test2 = mean(mean(Segment(seg).Forward(postomegafr))) > 0.5; end
            
            % (3,4) check F/R continuity before & after
            if seg > 1 test3 = (Segment(seg).BoundForward(1) == Segment(seg-1).BoundForward(2)); else test3 = NaN; end             
            if seg < numsegments test4 = (Segment(seg).BoundForward(2) == Segment(seg+1).BoundForward(1)); else test4 = NaN; end
            
            if nanmean([test1 test2 test3 test4])==0 & (max(fblock) < maxrevfr)
                Segment(seg).Forward = ~Segment(seg).Forward;
                Segment(seg).BoundForward = [Segment(seg).Forward(1), Segment(seg).Forward(Segment(seg).Length)];
                Segment(seg).Flip = 2; 
            end
            
            % Correct if there are too many reverse frames at the beginning, and a flip was added due to scrunch
            if any(scrunchflipfr)
                [ar,br] = findsegment(find(~Segment(seg).Forward));
                badflip = ar(find((br >= 30 | (br >= 10 & ar == br)) & ar <= max(scrunchflipfr)));
                if length(badflip) > 0
                    cr = addortho(scrunchflipfr,-badflip'); cr(find(cr < 1))=NaN;
                    undoflipfr = scrunchflipfr(find(sum(cr == repmat(nanmin(cr')',1,size(cr,2)))));
                    if length(undoflipfr) > 0
                        Segment(seg).Forward(1:undoflipfr(1)-1) = ~Segment(seg).Forward(1:undoflipfr(1)-1);
                        Segment(seg).BoundForward = [Segment(seg).Forward(1), Segment(seg).Forward(Segment(seg).Length)];
                        Segment(seg).Flip = 3;
                    end
                end
            end
            
            %if beginning of segment is no longer continuous with pervious segment, flip up to first omega
            if seg > 1 && (Segment(seg).BoundForward(1) ~= Segment(seg-1).BoundForward(2) & ...
                    abs(Segment(seg).PathAngVel(1)) < (Settings.MinRevAngVel / Settings.FrameRate))
                firstomegafr = min(find(Segment(seg).Omega | Segment(seg).Scrunch));
                if length(firstomegafr) > 0
                    Segment(seg).Forward(1:firstomegafr(1)-1) = ~Segment(seg).Forward(1:firstomegafr(1)-1);
                    Segment(seg).BoundForward = [Segment(seg).Forward(1), Segment(seg).Forward(Segment(seg).Length)];
                    if isfield(Segment,'Flip')
                        Segment(seg).Flip = Segment(seg).Flip + 0.1;
                    else
                        Segment(seg).Flip = 0.1;
                    end
                end
            end
                        
            % find number of omegas
            fwdhist = Segment(seg).Forward; 
            if seg > 1
                fwdhist = [Segment(seg-1).Forward, fwdhist]; 
                offset = length(Segment(seg-1).Forward);
            else
                offset = 0;
            end
            cleanfwd = imopen(imclose(fwdhist,se2),se2); % filter: fill gaps of 1-2; remove single or double fwd/rev
            [fend,flen] = findsegment(find(cleanfwd));
            [rend,rlen] = findsegment(find(~cleanfwd));
            Segment(seg).FwdBlock = [fend-flen+1; fend]' - offset; fend = fend - offset;
            Segment(seg).RevBlock = [rend-rlen+1; rend]' - offset; rend = rend - offset;

            [oend,olen] = findsegment(find(Segment(seg).Omega | Segment(seg).Scrunch));
            Segment(seg).OmegaBlock = [oend-olen+1; oend]';
            obegin = oend-olen+1; fbegin = fend-flen+1;
                        
            if (length(rend) == 0 & ((fbegin < 0 & fbegin > -5) | seg == 1)) rend = fbegin - 1; end
             
            realomegafr = [];
            if length(rend) > 0 & length(fend) > 0;
                fwdused = zeros(size(fend));
                for ob = 1:size(Segment(seg).OmegaBlock,1)
                    priorrev = (rend - oend(ob) < 0);
                    currentfwd = (fend >= obegin(ob)) & (fbegin-3 <= oend(ob)); % allow omega to occur 3 fr prior to fwd
                    
                    omegafr = any(priorrev) & any(currentfwd) & any(currentfwd ~= fwdused);
                    if omegafr
                        realomegafr = [realomegafr, obegin(ob)];
                        fwdused = fwdused | currentfwd; 
                    end
                end
            end
                        
            omegatooclose = find(diff(realomegafr) <= Settings.MinOmegaInterval*Settings.FrameRate);            
            realomegafr = realomegafr(setdiff(1:length(realomegafr),omegatooclose));
            
            Segment(seg).NumOmegas = length(realomegafr);
            Segment(seg).OmegaFr = Segment(seg).RealFrames(realomegafr);
            
            if Segment(seg).NumOmegas < 1
                Segment(seg).Pirouette = 0;
            end

        end
        
        [ae,be] = findsegment(find(Segment(seg).Erratic)); 
        if length(be)>0 ErraticLenMean = mean(be); else ErraticLenMean = 0; end
        if Segment(seg).Pirouette
            Segment(seg).ID = '7: Pirouette';
        elseif ErraticLenMean > (Settings.MinErrTime * Settings.FrameRate) %(sum(Segment(seg).Erratic) > (Settings.MinErrTime * Settings.FrameRate))
            Segment(seg).ID = '8: Erratic movement or swim';
        else
            Segment(seg).ID = '6: Long Turn or Reversal';
        end
    end
    
    Segment(seg).Code = str2num(Segment(seg).ID(1));
end

% forward turns
pathangle = [Segment.PathAngle];
pathangle2 = mod(pathangle + 180*~circshift([Segment.Forward]',1)',360);
dirangle = mod([Segment.PathAngle]+45,180)-45; % + 180*~circshift([Segment.Forward]',1)',360); % direction of fwd movement
dirangle2 = unwrap(dirangle * pi / 90) * 90/pi;
smdirang = smooth(dirangle2,7);
smdirang15 = round((smdirang-15)/60)*60 + 15;
cw = find(smdirang15(2:end) < smdirang15(1:end-1));
ccw = find(smdirang15(2:end) > smdirang15(1:end-1));

% Remove forward turns that are too close together (<= 4 frames?)
turnint = abs(addortho(cw,-ccw'));
[cwi,ccwi] = find(turnint <= 4);
validcw = []; validccw = [];
if length(cw)>0  validcw  = ~iselement(1:length(cw),cwi)'   & [1; diff(cw) > 4];  end
if length(ccw)>0 validccw = ~iselement(1:length(ccw),ccwi)' & [1; diff(ccw) > 4]; end

% Remove forward turns occurring in pirouettes or reversals
code = zeros(length(dirangle),1); for s=1:length(Segment), code(Segment(s).Frames) = Segment(s).Code; end
fwd = [Segment.Forward]; rev = ~fwd;
if length(cw)>0  validcw = validcw & iselement(code(cw),1:6) & fwd(cw)'; end
if length(ccw)>0 validccw = validccw & iselement(code(ccw),1:6) & fwd(ccw)'; end

%---------------
omegafr = []; % nrevfr = []; priorrev = [];
if isfield(Segment,'OmegaFr') omegafr = [Segment.OmegaFr]; end

reversal(nostallfr) = rev; reversal(stallfr) = 0;

SegAnalysis.Omega = [omegafr', PriorReversal(reversal,omegafr', ...
    (-(Settings.ReversalHistory * Settings.FrameRate):-1))];

SegAnalysis.FTurnCCW = [nostallfr(ccw(validccw))', PriorReversal(reversal,nostallfr(ccw(validccw))', ...
    (-(Settings.ReversalHistory * Settings.FrameRate):-1))];
SegAnalysis.FTurnCW = [nostallfr(cw(validcw))', PriorReversal(reversal,nostallfr(cw(validcw))', ...
    (-(Settings.ReversalHistory * Settings.FrameRate):-1))];

SegAnalysis.StallFr = stallfr';
SegAnalysis.NoStallFr = nostallfr;
SegAnalysis.CollFr = find(collision);
SegAnalysis.Length = length(CurrentTrack.Frames);
SegAnalysis.MaxFr = max(CurrentTrack.Frames);

SegAnalysis.OriginalDistance = origdist;

%--------------------------------------------------------------------------

function output = PriorReversal(reversals, framelist, frange) 

output = [];
if length(framelist)>0
    priorfr = addortho(framelist,frange);
    priorrev = double(reversals(max(1,priorfr))); 
    priorrev(find(priorfr < 1)) = NaN;
    
    vfr = sum(~isnan(priorrev'));
    nrevfr = nansum(priorrev') ./ (vfr + eps) * length(frange);
    output = [nrevfr',vfr'];
end

%--------------------------------------------------------------------------

function [ang,angvel,angacc,dist] = angle(pos)

len = size(pos,1);
if len >= 2
    smpos = [pos(1,:); (pos(2:len,:)+pos(1:len-1,:))/2; pos(len,:)];

    delpos = smpos(2:len+1,:)-smpos(1:len,:);
    ang = -atan2(delpos(:,2),delpos(:,1)) * 180 / pi;

    angvel = [0; ang(2:len)-ang(1:len-1)];
    angvel = mod(angvel+180,360)-180;

    angacc = [0; angvel(2:len)-angvel(1:len-1)];

    dist = sqrt(sum(delpos.^2,2));
else
    ang = 0;
    angvel = 0;
    angacc = 0;
    dist = 0;
end
%--------------------------------------------------------------------------

function [ang,angvel,angacc,dist] = angle2(pos)

len = size(pos,1);
if len >= 2
    delpos = pos(2:len,:)-pos(1:len-1,:);
    ang = -atan2(delpos(:,2),delpos(:,1)) * 180 / pi;
    ang = [ang(1); ang];

    angvel = [ang(2:len)-ang(1:len-1); 0];
    angvel = mod(angvel+180,360)-180;

    angacc = [0; angvel(2:len)-angvel(1:len-1)];

    dist = sqrt(sum(delpos.^2,2));
   dist = [0; dist];
else
    ang = 0;
    angvel = 0;
    angacc = 0;
    dist = 0;
end