function    makeLetterXColorHist(labelmatrix)

% function    makeLetterXColorHists()
% label matrix is an nx26 set of labels
% labels range from 0-11
letters = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

names = {...
    'black' ...
    'white' ...
    'red' ...
    'green' ...
    'yellow' ...
    'blue' ...
    'brown' ...
    'purple' ...
    'pink' ...
    'orange' ...
    'grey'...
    'none',...
    };

histcolors = [ 0 0 0;
    1 1 1;
    1 0 0;
    0 1 0;
    1 1 0;
    0 0 1;
    .8 .4 .12;
    .8 0 .8;
    1 .1 .55;
    1 .6 0;
    .5 .5 .5;
    0 0 0;
    ];




% let's do it again as percent;  adjusted for no matches
figure('name','percent of times each letter is given a color label', 'Color', [1 1 1],'Position',get(0,'ScreenSize'));
for i = 1:26
    
    subplot(6,5,i);
    [counts, bins]=hist(labelmatrix(:,i),12);
    counts  = counts/length(labelmatrix);
    hBar =bar(bins, counts,'hist');
    ylabel('% of subs');
%     set(gca,'XLim',[0 11],'XTick',bins,'XTickLabel',names,...
%         'YLim',[0 .8],'YTick',[0:.2:.6],'FontSize',6,'FontWeight','bold');
    set(gca,'XLim',[0 11],'XTick',[],'YLim',[0 .6],'YTick',[0:.2:.6],'FontSize',6,'FontWeight','bold');
    
    set(hBar,'FaceVertexCData',histcolors);
    title(letters(i),'FontSize',15,'FontWeight','bold');
    box off;
    
end




names = {...
    'black' ...
    'white' ...
    'red' ...
    'green' ...
    'yellow' ...
    'blue' ...
    'brown' ...
    'purple' ...
    'pink' ...
    'orange' ...
    'grey'...
    'none',...
    };

% need sorted vector for colors
% bwgrobrygblppn

huesorted = [ 1 11 2 3 10 7 5 4 6 8 9 12];

% sort letters by frequency

fqorder = [5 20  1 15  9 14 19  8 18  4 12 3 21  13 23 6 7 25  16  2 22   11 10 24 17 26];
% 
% let's do it again as percent;  adjusted for no matches
figure('name','percent of times each letter is given a color label sort colors', 'Color', [1 1 1],'Position',get(0,'ScreenSize'));
for i = [5 20  1 15  9 14 19  8 18  4 12 3 21  13 23 6 7 25  16  2 22   11 10 24 17 26];

    
    subplot(6,5,i);
    [counts, bins]=hist(labelmatrix(:,fqorder(i)),12);
    counts  = counts/length(labelmatrix);
    hBar =bar(bins, counts(huesorted),'hist');
    ylabel('% of subs');
%     set(gca,'XLim',[0 11],'XTick',bins,'XTickLabel',names,...
%         'YLim',[0 .8],'YTick',[0:.2:.6],'FontSize',6,'FontWeight','bold');
    set(gca,'XLim',[0 11],'XTick',[],'YLim',[0 .6],'YTick',[0:.2:.6],'FontSize',6,'FontWeight','bold');
    
    set(hBar,'FaceVertexCData',histcolors(huesorted,:));
    title(letters(fqorder(i)),...
        'FontSize',15,'FontWeight','bold');
    box off;
    
end







end
