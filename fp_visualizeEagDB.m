% fp_visualizeEagDB.m

% script which loads and then visualizes the eagleman database of
% grapheme-color synesthetes



% place to save figures
datavisdir = 'databasevisualizations';

% if dir doesn't exist, make it
if ~exist(datavisdir,'dir')
   mkdir(datavisdir); 
end



% load data
% this seems like it is happening twice
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
% figure('name', 'all matches in eagleman database', 'Color', [1 1 1]);
% imagesc(p_rgb);
% ylabel('SUBJECTS');
% xlabel('LETTERS');
% set(gca,'XTick',[1:26],'XTickLabel',labels, 'TickDir','out', 'YDir','normal');
% box off;
% title('All letter-color matches in eagleman database');
%


% let's assign a random rgb value to all the nans in order to break up any
% structure that comes from having nans all be black.

% get index to nans

indx =find(isnan(p_rgb));


% make a random rgb matrix of length index ,3
randcolors = rand(length(indx),1);

randp_rgb=p_rgb;
% fill in the nans with random colors
randp_rgb(indx) =randcolors;



% now we can look at the whole data set
figure('name', 'all matches in eagleman database', 'Color', [1 1 1],'Position',get(0,'ScreenSize'));
subplot(1,2,1);
imagesc(p_rgb);
ylabel('SUBJECTS');
xlabel('LETTERS');
set(gca,'XTick',[1:26],'XTickLabel',labels, 'TickDir','out', 'YDir','normal');
box off;
title('nans are black');

subplot(1,2,2);
imagesc(randp_rgb);
ylabel('SUBJECTS');
xlabel('LETTERS');
set(gca,'XTick',[1:26],'XTickLabel',labels, 'TickDir','out', 'YDir','normal');
box off;
title('nans are a random color');

saveas(gcf,[datavisdir '/AllMatchesAllSubs.png'],'png');
plot2svg([datavisdir '/AllMatchesAllsubs.svg'],gcf);

close(gcf);



% let's sort the letters by frequency to compare to the cognition paper
% under review from Herman et al.

% lewand frequency
% 5 20  1 15  9 14 19  8 18  4 12 3 21  13 23 6 7 25  16  2 22   11 10 24 17 26
% e t   a o   i n  s   h r   d l  c u   m  w  f g y p   b v    k  j  x  q  z
fqorder = [5 20  1 15  9 14 19  8 18  4 12 3 21  13 23 6 7 25  16  2 22   11 10 24 17 26];


% now we can look at the whole data set
figure('name', 'all matches in eagleman database sorted by frequency', 'Color', [1 1 1],'Position',get(0,'ScreenSize'));
subplot(1,2,1);
imagesc(p_rgb(:,fqorder,:));
ylabel('SUBJECTS');
xlabel('LETTERS');
set(gca,'XTick',[1:26],'XTickLabel',labels(fqorder), 'TickDir','out', 'YDir','normal');
box off;
title('fq in english');

% fq of first letter in a word

% 20 1 19 8 23 9 15 2 13 6 3 12 4 16 14 5 7 18 25 21 22 10 11 24 26 24
% t  a s  h w  i o  b m  f c l  d p  n  e g r  y  u  v  j  k  q  z  x

fqfrst = [20 1 19 8 23 9 15 2 13 6 3 12 4 16 14 5 7 18 25 21 22 10 11 24 26 24];

subplot(1,2,2);
imagesc(p_rgb(:,fqfrst,:));
ylabel('SUBJECTS');
xlabel('LETTERS');
set(gca,'XTick',[1:26],'XTickLabel',labels(fqfrst), 'TickDir','out', 'YDir','normal');
box off;
title('fq as first letter in a word in english');



% can see some cool stuff already.  for example there are numerous letters
% which seem to be strongly associated with colors overall.
%  A : RED
%  B : BLUE
%  C : YELLOW
%  D : DARK?
%  G : GREEN?
%  H : PINK/ORANGE?
%  I : WHITE
%  L : LIGHT?
%  O : WHITE
%  P : PINK/PURPLE?
%  R : RED
%  S : YELLOW?
%  X : BLACK
%  Y : YELLOW


% what we would like to see is the histograms across subjects for each
% letter.  one way would be to convert each rgb value to hsv and then do
% the histogram around the hue circle.  this is on the to do list.
% unfortunately it rules out the white letters

%
% % a second way is to plot the points on 3d plots where the axes are r,g,
% % and b, and look at if and how they cluster.
% % for each letter make a plot (later to be closed and saved as figures)
% i=11;
% % for i=1:1
% figure('name', ['matches for letter : ', labels(i)], 'Color', [1 1 1]);
% % what would be awesome would be to color each point by its rgb value
% clrs = squeeze(p_rgb(:,i,:));
% scatter3(p_rgb(:,i,1),p_rgb(:,i,2),p_rgb(:,i,3),20,clrs);
% xlabel('Red');ylabel('Green');zlabel('Blue');
% box off;
%
% % end
%
%
%
% % or we could plot every point in the data set just to see if it uniformly
% % covers the space (which it doesn't)
% figure('name','all colors generated in eagleman set across letters','Color',[1 1 1]);
% for i=1:26
%     clrs = squeeze(p_rgb(:,i,:));
%     scatter3(p_rgb(:,i,1),p_rgb(:,i,2),p_rgb(:,i,3),20,clrs);
%     hold on;
% end
% xlabel('Red');ylabel('Green');zlabel('Blue');
% box off;
%


%  a third way is to convert each of the rgb values to a label and then compute
% make some histograms showing counts or proportions in each bin
% first we need to convert our rgbs to color labels.  At this point the
% conversion is done by lookup table.  the lookup table was generated by
% manually labelling thousands of colors (a grid that overlays the rgb
% space).  we interpolate to the nearest point to get the color.  at the
% moment this lookup seems to be 90-95% accurate so it could be better.
% that set of labels is called lRGBnathan.mat

% let's make our labeled database
% this function takes the rgb colors and converts them to numerical labels
% and converted to color names
[dbNumbered dbLabeled] = fpRGB2ColorsJW(rgb);


%   dbLabeled            6588x26                11947850  cell
%   dbNumbered           6588x26                 1370304  double

% get distribution of color choices for each letter
makeLetterXColorHist(dbNumbered);


saveas(gcf,[datavisdir '/LetterxColorHistsAllSubs.png'],'png');
plot2svg([datavisdir '/LetterxColorHistsAllSubs.svg'],gcf);

close(gcf);

%
% fsize=get(gcf,'Position');
% set(gcf,'Position',[5 5 5*fsize(3) 5*fsize(4)]);
% saveas(gcf,'svgplots/letterHistsAll.jpg','jpg');
%




format short g;

% want to slip a table in here too in the command window
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(:,i),11);
    counts  = counts/length(dbLabeled);
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end

% would also like to convert the rgb to hsv.  matlab has a function for
% this, which I find absolutely never works correctly.  one problem is that
% the data contains some nans

% can see some cool stuff already.  for example there are numerous letters
% which seem to be strongly associated with colors overall.
%  A : RED
%  B : BLUE
%  C : YELLOW
%  D : DARK?
%  G : GREEN?
%  H : PINK/ORANGE?
%  I : WHITE
%  L : LIGHT?
%  O : WHITE
%  P : PINK/PURPLE?
%  R : RED
%  S : YELLOW?
%  X : BLACK
%  Y : YELLOW



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%   this section plots all the matches for each letter on 3d plots where
%   the axes are rgb.  as this is a lot of points these figures are slow to
%   make and adjust.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% % what we would like to see is the histograms across subjects for each
% % letter.  one way would be to convert each rgb value to hsv and then do
% % the histogram around the hue circle.  this is on the to do list.
% % unfortunately it rules out the white letters
%
%
% % a second way is to plot the points on 3d plots where the axes are r,g,
% % and b, and look at if and how they cluster.
% % for each letter make a plot (later to be closed and saved as figures)
% i=11;
% % for i=1:1
% figure('name', ['matches for letter : ', labels(i)], 'Color', [1 1 1]);
% % what would be awesome would be to color each point by its rgb value
% clrs = squeeze(p_rgb(:,i,:));
% scatter3(p_rgb(:,i,1),p_rgb(:,i,2),p_rgb(:,i,3),20,clrs);
% xlabel('Red');ylabel('Green');zlabel('Blue');
% box off;
%
% % end
%
%
%
% % or we could plot every point in the data set just to see if it uniformly
% % covers the space (which it doesn't)
% figure('name','all colors generated in eagleman set across letters','Color',[1 1 1]);
% for i=1:26
%     clrs = squeeze(p_rgb(:,i,:));
%     scatter3(p_rgb(:,i,1),p_rgb(:,i,2),p_rgb(:,i,3),20,clrs);
%     hold on;
% end
% xlabel('Red');ylabel('Green');zlabel('Blue');
% box off;

% these are too big to look at so maybe it would be helpful to smooth/bin the
% data

% counts of letters in each rounded 3d bin across all letters
allletterhists = [];

% rounded rgb values across all letters
allroundedclrs = [];
clrctr = 1;


% want to smooth the space a little bit
% one approach is simply to round towards some set of values
%     roundvec = 0:.25:1; %(1331 locations)
%     roundvec = 0:1/3:1;
%     roundvec = 0:.5:1;
    
% this introduces a problem because of edge effects.  so for example if your 
% points are 0, .5, and 1, then 0 gets stuff from the range 0-.25, .5 gets
% % .25-.75, and 1 gets .75-1.   since these are not absolute, you could
% normalize for the area on the fly


% 
% 
% % a second possibility is to do things histogram style and separate the
% % interval into equal bins and round to the center of those bins.  here you
% % don't have to normalize for the area, but since a lot of our data is at
% % the edges it may not show that off as well
% % here we would use linspace
% 
% numbins = 5;
% 
% % add 2 to numbins to get edges
% roundvec = linspace(0,1,numbins+2);
% 
% % now drop edges
% roundvec = roundvec(2:end-1);
% 
% 
% 
% % for each letter
% for i=1:26
%     clrs = squeeze(p_rgb(:,i,:));
%     %  all rgb values are between 0 and 1
%     % we could round to the nearest .1
%     % use downloaded function roundtowardvec
%     roundedclrs = roundtowardvec(clrs,roundvec);
% 
%     
%     allroundedclrs = [allroundedclrs; roundedclrs];
%     
%     % want to be able to scale our plots (either alpha or markersize) using the
%     % number of datapoints in each bin.  3d hist code in matlab looks like a
%     % pain so we could do it in a loop?
%     
%     counter=1;
%     bbplotmatrix = zeros(length(roundvec)^3,4);
%     
%     for a=1:length(roundvec)
%         for b=1:length(roundvec)
%             for c=1:length(roundvec)
%                 bbplotmatrix(counter,:) = [roundvec(a) roundvec(b) roundvec(c) ...
%                     sum(ismember(roundedclrs, [roundvec(a) roundvec(b) roundvec(c)],...
%                     'rows'))];
%                 counter=counter+1;
%             end
%         end
%     end
%     % less than a minute
% 
%     
% %     can add back in letter plots when we need them
%     
% %     figure('name',['distribution across rgb space for ' letters{i}],'Color',[1 1 1]);
% %     
% %  
% %     
% %     sizes=bbplotmatrix(:,4);
% %  
% %     
% %     % need to get rid of zeros
% %     p = find(sizes~=0);
% %     
% %     %     BUBBLEPLOT3(x,y,z,r,c), where c is a rgb-triplet array (in [0,1])
% %     %     with numel(x) rows, plots bubbles with colours specified by c.
% %     
% %     BUBBLEPLOT3(bbplotmatrix(p,1)',bbplotmatrix(p,2)',bbplotmatrix(p,3)',(sizes(p)/sum(sizes))',bbplotmatrix(p,1:3));
% %     % scatter3(roundedclrs(:,1),roundedclrs(:,i,2),roundedclrs(:,i,3))
% %     
% %     
% %     xlabel('Red');ylabel('Green');zlabel('Blue');
%     
%     %     collect everything into a giant matrix
%     
%     allletterhists = [allletterhists; bbplotmatrix];
%     
%     
% end
% 




% looks odd, because there are almost never strong values for blue
% have to sum across all the letters to get the whole distribution

% % make bubble plot for everything
% 
% figure('name',['distribution across rgb space for all letters'],'Color',[1 1 1]);
% 
% %    SCATTER3(X,Y,Z,S,C) displays colored circles at the locations
% %     specified by the vectors X,Y,Z (which must all be the same size).  The
% %     area of each marker is determined by the values in the vector S (in
% %     points^2) and the colors of each marker are based on the values in C.  S
% %     can be a scalar, in which case all the markers are drawn the same
% %     size, or a vector the same length as X,Y, and Z
% 
% sizes=allletterhists(:,4);
% 
% % will need to scale size
% % try % of max
% % s = bplotmatrix(:,4)/max(bplotmatrix(:,4));
% 
% % need to get rid of zeros
% p = find(sizes~=0);
% 
% %     BUBBLEPLOT3(x,y,z,r,c), where c is a rgb-triplet array (in [0,1])
% %     with numel(x) rows, plots bubbles with colours specified by c.
% 
% BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',10*(sizes(p)/sum(sizes))',allletterhists(p,1:3));
% 
% % values range from 1 to 1898
% % try using log scaling
% % this works but hard to see pattern
% % BUBBLEPLOT3(allletterhists(p,1)',allletterhists(p,2)',allletterhists(p,3)',(log(s(p))/(10*max(log(s(p)))))',allletterhists(p,1:3));
% 
% 
% % scatter3(roundedclrs(:,1),roundedclrs(:,i,2),roundedclrs(:,i,3))
% 
% 
% xlabel('Red');ylabel('Green');zlabel('Blue');
% 
% 
% 
% 
% 
% 
% % pair it with a uniform distribution
% figure('name',['uniform distribution across rgb space'],'Color',[1 1 1]);
% 
% % total number of matches
% numletterswithmatches = sum(sizes(p));
% % number for each point
% u = numletterswithmatches/(length(roundvec)^3);
% 
%  bplotmatrix(:,4) = u;
% 
%  BUBBLEPLOT3(bplotmatrix(:,1)',bplotmatrix(:,2)',bplotmatrix(:,3)',(10*bplotmatrix(:,4)/numletterswithmatches)',bplotmatrix(:,1:3));
% xlabel('Red');ylabel('Green');zlabel('Blue');
% 

% 
% % let's see how all letters are distributed in rgb
% % find all non nan rows
% p=find(~isnan(allroundedclrs(:,1)));
% 
% 
% figure('Name','histograms of rgb channels across letters','Color',[1 1 1]);
% subplot(3,1,1);
% [r bins] =hist(allroundedclrs(p,1),roundvec);
% bar(bins, r/sum(r), 'barwidth',1);
% box off;
% xlabel('red');
% ylabel('number of letters');
% subplot(3,1,2);
% [g, bins]=hist(allroundedclrs(p,2),roundvec);
% bar(bins, g/sum(g), 'barwidth',1);
% box off;
% xlabel('green');
% ylabel('number of letters');
% subplot(3,1,3);
% [b, bins] = hist(allroundedclrs(p,3),roundvec);
% bar(bins, b/sum(b), 'barwidth',1);
% box off;
% xlabel('blue');
% ylabel('number of letters');


% let's see how all CHROMATIC matches are distrbuted in rgb
% remove rows that have equal values in all three channels
% have already removed nanrows
% 
% chromrows = find(sum(allroundedclrs(p,:),2) ~= 3*allroundedclrs(p,1));
% 
% 
% 
% figure('Name','histograms of rgb channels across chromatic matches','Color',[1 1 1]);
% subplot(3,1,1);
% [r bins] =hist(allroundedclrs(chromrows,1),roundvec);
% bar(bins, r/sum(r), 'barwidth',1);
% box off;
% xlabel('red');
% ylabel('% of matches');
% subplot(3,1,2);
% [g, bins]=hist(allroundedclrs(chromrows,2),roundvec);
% bar(bins, g/sum(g), 'barwidth',1);
% box off;
% xlabel('green');
% ylabel('% of matches');
% subplot(3,1,3);
% [b, bins] = hist(allroundedclrs(chromrows,3),roundvec);
% bar(bins, b/sum(b), 'barwidth',1);
% box off;
% xlabel('blue');
% ylabel('% of matches');
% 
