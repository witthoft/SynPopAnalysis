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
    set(gca,'XLim',[0 11],'XTick',bins,'XTickLabel',names,...
        'YLim',[0 .8],'YTick',[0:.2:.6],'FontSize',6,'FontWeight','bold');
    set(hBar,'FaceVertexCData',histcolors);
    title(letters(i),'FontSize',15,'FontWeight','bold');
    box off;
    
end




end