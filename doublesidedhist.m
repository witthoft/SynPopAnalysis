% want to make a histogram where there is one histogram facing upward and
% another facing downward

% make some data
uph = [0 0 0 1 1 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 4 4 4 5 6 7];
downh = [ 0  1 1 2 3 3 3 4 5 5 5 5 5 5 5 5 6 6 6 6 6 7 7 7 7];


% my only idea is to make one histogram.

figure('Color',[1 1 1]); 

% we will need to control the axes position

arect = get(gca,'position');

% then have an axis handle to plot to as the order will matter for making
% nice plots
m = get(gcf,'currentAxes');
% left bottom, width height
% to have this work with the subplot function we will need to have it be
% generic to whatever the current axes are
% we want to move the bottom halfway up and use half the height
newrect = [ arect(1) arect(2)+(arect(4)/2) arect(3) arect(4)/2];

% move the lower edge of the axis so that it takes up less than half the
% figure window
set(m,'position',newrect);
% then create a new set of axes, the same size as the current ones but
% which are immediately below it.
h=axes('position',[arect(1), newrect(2)-newrect(4) arect(3), newrect(4)]);
% then invert the yaxis in the new axes and plot the second histogram

set(h,'Ydir','reverse');
hold on;
hist(downh,7);
% turn off the x values

% then plot the top half
set(gcf,'currentAxes',m);
hist(uph,7);


box off;


% don't want x ticks
set(h,'XTick',[]);
set(m,'XTick',[]);


% set the y limit of both axes to be the bigger of the two
ymax = max([get(h,'YLim') get(m,'YLim')]);
set(m,'YLim',[0 ymax]);
set(h,'YLim',[0 ymax]);

set(h,'Xcolor',[1 1 1]);










