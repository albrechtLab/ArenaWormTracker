%
% Display binned behavioral data in time and space.
%
% USAGE:
%   DensityPlot(data, filename)
% 
%   data: data variable from WormDensity script, i.e. All, Track, or Beh(i)
%   filename: name to display on figure

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 12-Apr-2011 15:35:21 
%---------------------------- 
function DensityPlot(data, filename)
if ~exist('filename') || isempty(filename) filename = ''; end

densfig = findobj(get(0,'Children'),'Tag','Worm Density');
if isempty(densfig) densfig = figure; set(densfig,'Tag','Worm Density'); end
figure(densfig); clf;
dirlabel = {'up','horiz','down'};

hy = subplot(5,4,1:3); imagesc(data.Time,data.Y,data.TimeY'); colormap([flipud(gray(64));zeros(64,3)]); ylabel('y(mm)'); 
title(filename,'Interpreter','none');
hx = subplot(5,4,5:7); imagesc(data.Time,data.X,data.TimeX'); colormap([flipud(gray(64));zeros(64,3)]); ylabel('x(mm)'); xlabel('Time (min)');
subplot(5,4,4); barh(data.Y,sum(data.XY,1)); ylim(get(hy,'YLim')); axis ij; 
if isfield(data,'Text')
    title(data.Text);
else
    title('(Tracked animals)');
end
subplot(5,4,8); barh(data.X,sum(data.XY,2)); ylim(get(hx,'YLim')); axis ij;
n = 3; for j=1:3; for i=1:n; 
    subplot(5,n+1,8+((j-1)*(n+1))+i); 
    imagesc(data.X,data.Y,squeeze(msum(data.XYTimeD(:,:,1+round(end*(i-1)/n):round(i*end/n),j),3:4))'); 
    colormap([flipud(gray(64));zeros(64,3)]); axis equal tight; title(sprintf('%d/%d %s',i,n,char(dirlabel(j)))); if i==1 ylabel('y(mm)'); end; if j==3 xlabel('x(mm)'); end;
    end;
    subplot(5,n+1,8+j*(n+1)); imagesc(data.X,data.Y,squeeze(msum(data.XYTimeD(:,:,:,j),3:4))'); axis equal tight; title(sprintf('All %s',char(dirlabel(j)))); if j==3 xlabel('x(mm)'); end;
end; 

orient tall
