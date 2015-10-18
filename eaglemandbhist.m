% want to make a histogram where there is one histogram facing upward and


% another facing downward

% histogram
% colors come from fp_visualizeEagDB
figure('name','histogram of color names for each letter for no group syns', 'Color', [1 1 1]);
for i = 1:26
    
    subplot(5,6,i);
    %     get our data
    % fq
% there is something very retarded about the way bin widths are determined
% for example, suppose you have categories 0 to 10 which you want divided
% into 11 bins.  you do not get 11 bins of width 1 because matlab does not
% realize your data is categorical.  instead you get bins of widt 10/11.
% not sure what the reasoning is
% moreover, if you have categories 0:10 and two of those categories have 0
% counts in them, matlab makes the bins have widths of 9/11 and so even
% though you are just making a new histogram of the same kind of data, your
% bar widths are different.  so fuck their bins
    [counts, bins]=hist(dbNumbered(find(syntype==1),i),11)
    counts  = counts/length(find(syntype==1));
    %    m
%     let's do it for every synesthete who is not in the magnet set
%    [fcounts, fbins]=hist(dbNumbered(find(syntype~=1),i),11)
%     fcounts  = fcounts/length(find(syntype~=1));
% 
    [fcounts, fbins]=hist(dbNumbered(find(syntype==2),i),11)
    fcounts  = fcounts/length(find(syntype==2));
%     
    
    
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
    mBar =bar(.5:1:10.5, counts,'hist');
    set(mBar,'FaceVertexCData',histcolors);
    % turn off the x values
    % then plot the top half
    set(gcf,'currentAxes',m);
    
    fBar =bar(.5:1:10.5, fcounts,'hist');
    set(fBar,'FaceVertexCData',histcolors);
    box off;
    
    
    % don't want x ticks
    set(h,'XTick',[]);
    set(m,'XTick',[]);
    
    
    % set the y limit of both axes to be the bigger of the two
    %     ymax = max([get(h,'YLim') get(m,'YLim')]);
    %     set(m,'YLim',[0 ymax]);
    %     set(h,'YLim',[0 ymax]);
    
    % cute but for this we need it to be constant
    
    ymax = .5;
    set(m,'YLim',[0 ymax]);
    set(h,'YLim',[0 ymax]);
    
    set(h,'Xcolor',[1 1 1]);
    
    ylabel('% of matches');
    title(letters(i));
    
    %     need y = 0 linw
    hold on;
    set(gcf,'currentAxes',m);
    plot(0:.1:10,0,'-k','Linewidth',4);
    box off;
end








