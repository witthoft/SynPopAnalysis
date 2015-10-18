% fp_visualizeEagDB.m

% script which loads and then visualizes the eagleman database of
% grapheme-color synesthetes

% load data
load Engl_Alpha_5.mat;

% rgb has all our matching data
%   rgb                6588x3x26
% let's permute our rgb matrix so it is an image

p_rgb = permute(rgb, [1,3,2]);
%  p_rgb              6588x26x3


% unfortunately there are a few values which are just a tiny bit smaller
% than 0 and some which are a tiny bit larger than 1
% just reset those
p_rgb(p_rgb>1)=1;
p_rgb(p_rgb<0)=0;

% now we can look at the whole data set
figure('name', 'all matches in eagleman database', 'Color', [1 1 1]);
imagesc(p_rgb(:,nooverlap,:));
ylabel('SUBJECTS');
xlabel('LETTERS');
set(gca,'XTick',[1:length(nooverlap)],'XTickLabel',no_o_letters, 'TickDir','out', 'YDir','normal');
box off;
Title('All letter-color matches in eagleman database');


[dbNumbered dbLabeled] = fpRGB2ColorsJW(rgb);


%   dbLabeled            6588x26                11947850  cell
%   dbNumbered           6588x26                 1370304  double
% the lookup table for numbers to names is

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
    .5 .5 .5;];

% 
% 
% 
% % ok:  now we want some histograms
% % since we want to color each bar in the right color can't really use the
% % hist function
% figure('name','histogram of color names for each letter ', 'Color', [1 1 1]);
% for i = 1:26
%     
%     subplot(6,5,i);
%     [counts, bins]=hist(dbNumbered(:,i),11)
%     hBar =bar(bins, counts,'hist');
%     ylabel('Number of Subjects');
%     %    set(gca,'XTick',[0:length(names)],'XTickLabel',names,'YLim',[0
%     %    3000],'FontSize',6);
%     set(gca,'YLim',[0 3000],'FontWeight','bold');
%     set(hBar,'FaceVertexCData',histcolors);
%     title(labels(i));
%     box off;
%     
% end


% let's do it again as percent
figure('name','percent of times each letter is given a color label', 'Color', [1 1 1]);
for i = 1:26
    
    subplot(6,5,i);
    [counts, bins]=hist(dbNumbered(:,i),11)
    counts  = counts/length(dbLabeled);
    hBar =bar(bins, counts,'hist');
    ylabel('Number of Subjects');
    %    set(gca,'XTick',[0:length(names)],'XTickLabel',names,'YLim',[0
    %    3000],'FontSize',6);
    set(gca,'YLim',[0 .6],'FontWeight','bold');
    set(hBar,'FaceVertexCData',histcolors);
    title(labels(i));
    box off;
    
end



format short g;

% want to slip a table in here too in the command window
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(:,i),11);
    counts  = counts/length(dbLabeled);
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end



