%
% Bar plot by groups, showing means and std/err, with statistics 
%
% USAGE: output = barplotstat(datavector, groupvector, errbar, plist, pmarker, compgroup, comparisontype, ShowPoints)
%
%       datavector is vector of values
%       groupvector is vector of names assigning groups to each data value; same size as datavector
%       errbar is 'std' (default) or 'sem'
%       plist is list of p-values to test (default is : [0.05 0.01 0.001 0.0001]
%       pmarker is cell array of markers to label p-values on bars
%                                         (default: {'.','+','*','**'})
%       compgroup is name of a group for p-value comparison
%       comparisontype is 'tukey-kramer' (default), 'bonferroni', etc.
%       ShowPoints, if true, shows scatter of data points;
%           if vector, point styles are grouped by vector value

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 
% 11-Feb-2011 10:41:34 
%---------------------------- 

function output = barplotstat(datavector, groupvector, errbar, plist, pmarker, compgroup, comparisontype, ShowPoints)

if isnumeric(groupvector)
    if diff(size(groupvector))>0 groupvector = groupvector'; end
    groupvector = cellstr(num2str(groupvector));
end

if nargin < 8 || isempty(ShowPoints) ShowPoints = false; end
if nargin < 7 || isempty(comparisontype) comparisontype = 'tukey-kramer'; end
if nargin < 6 || isempty(compgroup) compgroup = unique(groupvector); end
if nargin < 5 || isempty(pmarker) pmarker = {'.','+','*','**'}; end
if nargin < 4 || isempty(plist) plist = [0.05 0.01 0.001 0.0001]; end
if nargin < 3 || isempty(errbar) errbar = 'std'; end

stat = anova1multicompare(datavector,groupvector,plist,comparisontype);

% calculate std, sem
for i = 1:length(stat.groups)
    index = find(strcmp(stat.groups(i),groupvector));
    stat.anovast.mean(i) = mean(datavector(index));
    stat.anovast.std(i) = std(datavector(index));
    stat.anovast.sem(i) = stat.anovast.std(i) / sqrt(stat.anovast.n(i));
end

if strcmp(errbar,'sem') 
    eb = stat.anovast.sem;
else
    eb = stat.anovast.std;
end

if isnumeric(compgroup) 
    compidx = compgroup; 
else
    [sortedgroups,sortidx] = sort(stat.groups);
    isortidx = []; for i=1:length(sortidx); isortidx(i)=find(sortidx==i); end
    compidx = find(strcmp(compgroup,sortedgroups));
end

h = barploterr(stat.anovast.means,eb,stat.groups);
ylim = get(gca,'YLim'); delta = 0.05*diff(ylim);

nbars = length(unique(groupvector));
ncol = max(5,nbars); sep = floor(254/(ncol-1));
cmap = jet(sep*(ncol) + 1); cmap = cmap(1:sep:size(cmap,1),:);
colormap(cmap(1:nbars,:));
barcolors = cmap;

% cmap = colormap;
% cmaps = size(cmap,1);
% barcolors = cmap(round(1:((cmaps-1)/(nbars-1)):cmaps),:)

pmarker = [{'ns'},pmarker];
ht = [];
for compn = 1:length(compidx)
    [m,n] = find(stat.stats(:,1:2) == compidx(compn));
    pdat = sortrows([diag(stat.stats(m,3-n)),m,stat.stats(m,4)]);

    if length(compidx)==1
        pmarkerpos = eb + stat.anovast.means + delta*compn;
        pmarkerpos(find(eb + stat.anovast.means < 0)) = delta*compn;
    else
        pmarkerpos = max(0,repmat(max(datavector),1,nbars)) + delta*compn;
    end

    ht = [ht; text(pdat(:,1),pmarkerpos(pdat(:,1)),char(pmarker(pdat(:,3)+1)), ...
          'HorizontalAlignment','center','Color',0.5*barcolors(compidx(compn),:))];
end

if ShowPoints == 1
    plot(grp2idx(groupvector),datavector,'ko','MarkerSize',4);
elseif length(ShowPoints)==length(datavector)
    [pa,pb,pc] = unique(ShowPoints);
    mlist = ['ko';'ks';'kd';'kx'];
    gi = grp2idx(groupvector);
    for i=1:max(pc)
        pi = find(pc == i);
        plot(gi(pi),datavector(pi),mlist(i,:),'MarkerSize',5);
    end
end

stat.pdat = pdat;
set(gca,'UserData',stat);

output = [h, ht'];



%--------------------------------------------------------------------------

function output = anova1multicompare(datavector, groupvector, plist, comparisontype)
% USAGE: output = anova1multicompare(datavector, groupvector, plist, comparisontype)
%
%       plist is list of p-values to test (default is : [0.05 0.01 0.001 0.0001]
%       comparisontype is 'tukey-kramer' (default), 'bonferroni', etc.
%

if nargin < 4 comparisontype = 'tukey-kramer'; end
if nargin < 3 plist = [0.05 0.01 0.001 0.0001]; end

plist = [1,sort(plist,'descend')];

ds = size(datavector);
gs = size(groupvector);
if ~all(ds == gs)
    if all(ds == fliplr(gs))
        gs = gs';
    else
        error('datavector and groupvector must have same size');
        return
    end
end
    
[p,t,st] = anova1(datavector,groupvector,'off');

output.anovap = p;
output.anovat = t;
output.anovast = st;

n = 0; for i = 1:length(st.means)-1; n = n+i; end
pvali = ones(n,1);
for i = 2:length(plist)

    [c,m,h,gnames] = multcompare(st,'alpha',plist(i),'ctype',comparisontype,'display','off');

    ptest = ~xor(c(:,3)>0, c(:,5)>0);  % true if comparison p-value at least plist(i)

    pvali = pvali + ptest;

    output.multcomp(i-1).alpha = plist(i);
    output.multcomp(i-1).c = c;
end

pval = plist(pvali);
ns = find(pvali == 1);

output.stats = [c(:,[1 2 4]), pvali-1, pval'];
if size(c,1) > 1
    output.statcell = [gnames(c(:,1:2)), num2cell(c(:,4)), num2cell(pval')];
else
    output.statcell = [gnames(c(:,1:2))', num2cell(c(:,4)), num2cell(pval')];
end
output.statcell(ns,4) = {'n.s.'};
output.groups = gnames;
    