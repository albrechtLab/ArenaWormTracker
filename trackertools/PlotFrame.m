function PlotFrame(FigH, Mov, Tracks)

warning('off','Images:initSize:adjustingMag')

figure(FigH)
if isfield(Mov,'cdata');
    imshow(Mov.cdata, Mov.colormap);
else
    imagesc(Mov);
end
hold on;

if ~isempty(Tracks)
    ActiveTracks = find([Tracks.Active]);
else
    ActiveTracks = [];
end

for i = 1:length(ActiveTracks)
    figure(FigH)
    plot(Tracks(ActiveTracks(i)).Path(:,1), Tracks(ActiveTracks(i)).Path(:,2), 'r');
    plot(Tracks(ActiveTracks(i)).LastCoordinates(1), Tracks(ActiveTracks(i)).LastCoordinates(2), 'b+');
end

drawnow;
hold off;    % So not to see movie replay
