%
% Analyze and correct specific behaviors in worm tracks.
%
% USAGE:
%   TrackAnalysis = AnalyzeTrack(Segment,SegAnalysis,trackbox,Settings)
% 
%   Segment: structure containing data for each track segment
%   SegAnalysis: structure of post-segmentation analysis
%   trackbox: dimensions of valid arena tracking bounds
%   Settings: segmentation settings structure
%   TrackAnalysis: output analysis structure

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 29-Mar-2011 11:40:00 
%---------------------------- 

function TrackAnalysis = AnalyzeTrack(Segment,SegAnalysis,trackbox,Settings)

if nargin < 4
    SegmentationSettings;
end

FRP(SegAnalysis.StallFr) = 3;
FRP(SegAnalysis.NoStallFr) = 1 + ~[Segment.Forward];

rev = (FRP == 2);

% reconstruct x-y path and animal heading across pauses
xy = zeros(length(FRP),2);
for s = 1:length(Segment)
    xy(Segment(s).RealFrames,:) = Segment(s).Path;
end

[a,b] = findsegment(SegAnalysis.StallFr);
for i=1:length(a)
    xy(a(i)-b(i)+1:a(i),:) = repmat(xy(a(i)-b(i),:),b(i),1);
end

%-----------------
%  Analyze turns
%-----------------

om = []; curveCW = []; curveCCW = [];

if (isfield(SegAnalysis,'Omega') && size(SegAnalysis.Omega,1) > 0) om = SegAnalysis.Omega(:,1:2); end
if (isfield(SegAnalysis,'FTurnCCW') && size(SegAnalysis.FTurnCCW,1) > 0) curveCCW = SegAnalysis.FTurnCCW(:,1:2); end
if (isfield(SegAnalysis,'FTurnCW') && size(SegAnalysis.FTurnCW,1) > 0) curveCW = SegAnalysis.FTurnCW(:,1:2); end

TrackAnalysis.TurnOmega = [ [curveCW, 1*ones(size(curveCW,1),1)]; [curveCCW, 3*ones(size(curveCCW,1),1)]; ...
                            [om, 5*ones(size(om,1),1)] ];

%-----------------------------------
%  Analyze frame-by-frame behavior
%-----------------------------------

ncode = []; nrev = [];
for seg = 1:length(Segment)
    ncode(SegAnalysis.NoStallFr(Segment(seg).Frames)) = Segment(seg).Code;
    nrev(SegAnalysis.NoStallFr(Segment(seg).Frames)) = ~Segment(seg).Forward;
end

[a,b] = findsegment(SegAnalysis.CollFr);
for ci = 1:size(a,1)
    ncode(a(ci)-b(ci)+1 : min(a(ci) + Settings.CollEffectTime * Settings.FrameRate, SegAnalysis.Length)) = 10;
end

oobfr = [];
if exist('trackbox')
    oobfr = xy(:,1) < (trackbox(1) + Settings.XBorderPadding) | xy(:,1) > (sum(trackbox([1 3])) - Settings.XBorderPadding);  % frames out of bounds
    oobfr = imdilate(oobfr',strel('line',(Settings.XBorderEffectTime*Settings.FrameRate)+1,0))';
end

% Correct short rev --> omega rev
nrevp = imclose(nrev(SegAnalysis.NoStallFr),strel('line',2,0)); % no stalls...
[a,b] = findsegment(find(nrevp));
for p = 1:length(a)
    frsR2O = SegAnalysis.NoStallFr(min(a(p)-b(p)+1:a(p)+2,length(SegAnalysis.NoStallFr))); %frsR2O = frsR2O(find(frsR2O <= length(ncode)));
    if length(frsR2O) > 0 && any(ncode(frsR2O) == 7) ncode(frsR2O) = 7; end
end

% Correct omega fwd before rev --> loopy fwd
[a1,b1] = findsegment(find(ncode(SegAnalysis.NoStallFr)==7)); %omegaseg = find([Segment.Code]==7);
for o = 1:length(a1)
    frs = SegAnalysis.NoStallFr(a1(o)-b1(o)+1:a1(o)); 
    orevs = nrev(frs);
    [a,b] = findsegment(find(orevs)); c = find(b > 1,1); firstrev = a(c)-b(c)+1;
    if firstrev > 1
        ncode(frs(1:firstrev-1)) = 6; % change to loopy fwd
    end

    [a2,b2] = findsegment(find(~orevs)); c2 = find(b2 > 20); 
    for of = 1:length(c2)
        longforward = frs(a2(c2(of))-b2(c2(of))+1:a2(c2(of)));
        longFxy = xy(longforward,:);
        if length(longFxy > 0)
            longFdist = sqrt(sum((longFxy - repmat(longFxy(1,:),length(longforward),1)).^2,2));
            ncode(longforward(longFdist > 0.5*Settings.PixelSize)) = 6;
        end
    end
end

%reanalyze corrected pirouettes
omegafr = []; if isfield(Segment,'OmegaFr') omegafr = [Segment.OmegaFr]; end
[ap,bp] = findsegment(find(ncode(SegAnalysis.NoStallFr)==7));
if length(ap)>0
    pblock = SegAnalysis.NoStallFr([ap-bp+1; ap]');
    yespir = zeros(size(pblock,1),1);
    for pb = 1:size(pblock,1)
        if ~any(omegafr >= pblock(pb,1) & omegafr <= pblock(pb,2))
            ncode(pblock(pb,1):pblock(pb,2)) = 6; % change to loopy fwd
        end
    end        
end

ncode(SegAnalysis.StallFr) = 9;
ncode(oobfr) = 10;
nrev(SegAnalysis.StallFr) = 0;
nFRP = 1 + nrev; nFRP(SegAnalysis.StallFr) = 3;

TrackAnalysis.Code = ncode;
TrackAnalysis.FRP = nFRP;
TrackAnalysis.Beh = diag(Settings.BehDef(nFRP,ncode))';

if size(TrackAnalysis.TurnOmega,1) > 0
    TrackAnalysis.TurnOmega = [TrackAnalysis.TurnOmega, TrackAnalysis.Code(TrackAnalysis.TurnOmega(:,1))', ...
                               xy(TrackAnalysis.TurnOmega(:,1),:), ];
    TrackAnalysis.TurnOmega(:,3) = TrackAnalysis.TurnOmega(:,3) + (TrackAnalysis.TurnOmega(:,2) >=  Settings.LongReversalTime*Settings.FrameRate);
    TrackAnalysis.TurnOmega = TrackAnalysis.TurnOmega(find(TrackAnalysis.TurnOmega(:,4) ~= 10),:); % remove oob turns
end