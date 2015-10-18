
% wrapper for analyzing the eagleman database generally

% data and functions on my laptop
cd ~/Dropbox/synesthesia/Eagleman/

% this is the data

load EaglemanColoredAlphabets.mat


% whos
%   Name            Size                Bytes  Class     Attributes
%
%   READ_ME         1x106                 212  char
%   ans             1x49                   98  char
%   labels         26x1                    52  char
%   u_rab        6588x3x26            4110912  double
%   u_rgb        6588x4x26            5481216  double
%   u_rlab       6588x4x26            5481216  double
%   userid       6588x1                 52704  double


% here is a bunch of additional information provided as a function
% AboutPtCols'

% ans =
%
%     'users_id'
%     'BatteryId'
%     'txtAge'
%     'optGender'
%     'txtMotherTongue'
%     'optSynInFamily'
%     'optHandedness'
%     'optMusicPitch'
%     etc....

LoadAboutPt

% whos
%   Name                 Size                  Bytes  Class     Attributes
%
%   AboutPt          54450x31              107375498  cell
%   AboutPtCols          1x31                   2674  cell
%   this file is much bigger than the whole set as it is info on everyone
%   in the eagleman database.  somewhere there is code that matches these,
%   though I think to get the ages right I had to do some hand editing as
%   people entered the age in very inconsistent formats





% this is handy for labeling figures
letters = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};

% color names in the order we number them when we label all the colors in
% the database as categories
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


% useful for coloring figures that use color categories
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




% show me the nans in the database
s_nansdistribution




% visualize the whole database
% and % color category for each letter
fp_visualizeEagDB








% suppose we wanted to something weird, like estimate the expected
% proportions of colors from the labeled database.  that is, what
% percentage of the space is black etc....
load lRGBnathan.mat;

per_colorspace=[];
for i=0:11
    per_colorspace = [per_colorspace length(find(labeledRGB == i))];
end



% check that it adds up to 972
% sum(per_colorspace)
disp('labeled rgb space color proportions');
for i=1:length(names)
    disp([names{i} '    '  num2str(per_colorspace(i)/sum(per_colorspace),'%0.2f')]);
end



% then you would want to see the % of matches that are in each category


% %
% just brings out the things that can be seen from the full dataset.
% A is red more than 50% of the time
%  B is blue more than 40% of the time
% c is yellow around 38% of the time
%  D has no obvious hue preference
%  E does seem to be dominated by blues yellows and greens not sure why.
% G is green about 37% of the time
%  h is orange about 30% of the time
%  I is black or white about 58% of the time but mostly black
%  L has a slight tendency towards yellow or blue
%  M is red about 30% of the time
%  N is orange about 22 % of the time
%  O is black or white about 60 % of the time
%  R is red about 38% of the time
% X is blacke about 39% of the time
%  Y is yellow about 42% of the time
%    Z is black about 38% of the time




% next we would like to find out how many synesthetes in the database were
% likely to have had the magnet set.  running this script generates a
% number of useful figures

% s_fpEaglemanMatches

s_fpRichAndEaglemanMatches


% what proportion of all matches have each label

alllabels = reshape(labels.eagleman,...
    size(labels.eagleman,1)*size(labels.eagleman,2),1);

% remove unmatched letters
matchedletterlabels = alllabels(find(alllabels~=11));

% get %s
[countcatmatches catbins] = hist(matchedletterlabels,unique(matchedletterlabels));

pctcatmatches = countcatmatches/sum(countcatmatches);

% see how out of whack matches are with color space proportions
fprintf('label\tspace\tmatches\tratio\n');
for i=1:length(names)-1
    fprintf('%s\t%4.2f\t%4.2f\t%4.2f\n',...
        names{i},per_colorspace(i)/sum(per_colorspace),pctcatmatches(i),...
        (pctcatmatches(i))/(per_colorspace(i)/sum(per_colorspace)));
end


% make these into bar plots

figure('Name','proportions of color space','Color',[1 1 1 ]);
% proportion of rgb space with a given label
subplot(3,1,1);
hBar=bar(0:10,per_colorspace(1:end-1)/sum(per_colorspace),'hist');
box off
set(gca,'XTickLabel',names,'XLim',[-0.5,10.5],'YLim',[0,.5]);
ylabel('% of rgb space');
set(hBar,'FaceVertexCData',histcolors(1:end-1,:));

% proportion of responses with a given label
subplot(3,1,2);
hBar=bar(0:10,pctcatmatches,'hist');
box off
set(gca,'XTickLabel',names,'XLim',[-0.5,10.5],'YLim',[0,.5]);
ylabel('% of responses');
set(hBar,'FaceVertexCData',histcolors(1:end-1,:));


% why not just as scatterplot
figure('Name','% of RGB space vs. % of responses','Color',[1 1 1]);
% the scatterplot function sucks for coloring so use loop as it is short
pctspace = per_colorspace(1:end-1)/sum(per_colorspace);

% plot points in each color
for i=1:10
    plot(pctspace(i),pctcatmatches(i),'o','MarkerSize',20,'MarkerEdgeColor','k',...
        'MarkerFaceColor',histcolors(i,:));
    hold on
end

% plot x=y line
plot(0:0.01:.35,0:0.01:.35,'k--');

box off;
set(gca,'XLim',[0,.35],'YLim',[0,.35]);
axis square;
xlabel('% of RGB space');
ylabel('% of responses');




% could work out how overrepresented each point in our binned space is
% how about bars where they are number observed / uniform and ordered and
% colored.  1 would mean that there are exactly as many matches as expected
% given the uniform distribution  this should also be a function









% to obtain an actual probability, let's shuffle the data set many times
% and count the probability of observing n or more matches

s_pNorMoreMatchestoMagnets



%  the first three figures show variations on the number of subjects with n matches to the
%  letter set ranging from 0 to 26 letters.   this was done by assigning
%  each match to a color label by comparing it to a labeled rgb cube (which
%  might still need some work).  the green line shows the number of
%  subjects with n matches to the set when the data is shuffled. the two
%  lines begin to diverge by 8 letters but there are only 10 subjects with
%  10 matches to the letter set when the data is shuffled, but 100 when not
%  shuffled for the same number of matches.  there are no subjects which
%  have more than 10 matches when the shuffled data set is compared to the magnets
% while there are hundreds when using the unshuffled data

% there are about 400 or so if one strictly counts more than 10 matches.
% of course even 8 matches is unlikely, at least if one views the problem
% as a series of independent random draws with 11 choices and with
% replacement.  so the question is, is the green line a fair null
% distribution?  one question is if it matters that some of the bins
% repeat.  i.e. there are 4 red letters.  given that a sizeable percentage
% of the subjects appear to have had the magnet set, does that skew the
% green distribution rightward?  that might happen because those subjects
% have ideally only 6 different colors of letters and also only 6 bins.  in
% any case, the number of matches expected in that case would be much
% higher than when thought about as a binomial distribution.

% the fourth figure shows that the real data and the shuffled data show the
% same distribution of matches to a shuffled data set.  which is to be
% expected.

fp_visualizeShuffledEagDB

% or, suppose we pulled out those that have more than 10 matches to the
% fridge magnets and shuffled the remainder.  what does that distribution
% look like?  have to make up those 500 subjects unless the curves are
% normalized.


% another thing is to run this same analysis but comparing subjects to the
% strong group level frequencies in the data.

% s_gtEaglemanMatches
s_fqEaglemanMatches


% this also reminds me that perhaps we should look for subjects who give
% the same color for most letters.  that would be a likely sign of fraud or
% at least unverifiability, and would impact the shuffling in a weird way
% if there were a lot of them.


% we should probably see if we can assign each synesthete to a group and
% generate separate plots for each.  really the point is to analyze the
% subjects left after you take out the frequent template group and the
% magnet group

syntype = zeros(length(dbLabeled),1);


% num of matches needed to be in group
magnetthreshold = 10;
fqthreshold = 10;

% code for labels variable is in s_fpEaglemanMatches but in this script
% should be in workspacelabels.
magmatches = labels.eagleman == labels.magnet;
fqmatches = labels.eagleman == labels.fq;
% let's set the fq first
syntype(find(sum(fqmatches,2)>=fqthreshold))=1;
% then magnet.
syntype(find(sum(magmatches,2)>=magnetthreshold))=2;
% maybe want to find those that are in both groups?


% how overrepresented is each point in rgb 
f_colorspaceRepBars(p_rgb,roundvec,'allmatches')



% visualize matching distributions for just the magnet synesthetes (same as
% code for fp_visualizeEagDatabase so should be a function
% make bubble plot for everything


% define rounding
% in function above...fix

% pull out magnet subset
magnetmatches = p_rgb(find(syntype==2),:,:);

% make bubbleplots
f_colorspacebubbleplot(magnetmatches,roundvec,'magnetsyns');
f_colorspaceRepBars(magnetmatches,roundvec,'magnetsyns')


% pull out only nonmagnetsynesthes
nonmagnetmatches = p_rgb(find(syntype~=2),:,:);
f_colorspacebubbleplot(nonmagnetmatches,roundvec,'not magnetsyns');
f_colorspaceRepBars(nonmagnetmatches,roundvec,'not magnetsyns')




% %
% need to create a distance matrix for each subject
% will use euclidian to start with
% can do this for each subject using pdist




% for 26 letters there are 26^2-26/2=325 comparisons

% letter pairs in order of distance measures using order that is output by
% pdist
letterpairs = cell(325,1);
counter=1;
for i=1:26
    for j=1:26
        if j>i;
        letterpairs{counter} = [letters(i) letters(j)]; 
       counter = counter+1;
        end
    end
end
    



% matrix which will hold all our distances
% m comparisons by n subjects
rgb_dists = nan(325,6588);


% calculate distance matrix for each subject
for i=1:length(p_rgb)
%     fill with calculated pdists
% comparisons with nans appear as nans
    rgb_dists(:,i)=pdist(squeeze(p_rgb(i,:,:)));

end

% visualize average distance matrix

% meanrgbdist = nanmean(rgb_dists,2);
% try median
meanrgbdist = nanmedian(rgb_dists,2);

 figure('Name','histogram of average distance between pairs of letters in rgb','Color',[1 1 1]);
 
 hist(meanrgbdist)
 xlabel('euclidian distance in rgb');
 ylabel('count of letter pairs');
 
%  plot against linear distance in sequence

% 
% this will be handy for other things
letterseqdists = pdist((1:26)');


figure('Name','rgb distance vs alphabet distance','Color',[1 1 1]);

plot(letterseqdists,meanrgbdist,'ro');
box off;
xlabel('alphabetic distance');
ylabel('mean rgb distance');

% looks like it doesn't matter much for english
% add in correlation stuff later.


% here are the non gibson similarity measures kindly loaned to us by marcus
% watson
simmeasures = csvread('Non-GibsonLetterDistances.csv',1, 2);

gibsonshapemeasure = [2.645751311
2.449489743
2.645751311
2.449489743
2.236067977
2.449489743
1.732050808
2
2.828427125
1.414213562
2.645751311
2
2.236067977
2.449489743
2.449489743
2.236067977
2.236067977
2.828427125
2
2.449489743
1.732050808
2
1.414213562
1.732050808
2.236067977
2.236067977
1.414213562
1.732050808
2.449489743
2.645751311
2
1.732050808
2.645751311
2.236067977
2.449489743
2.236067977
2.449489743
1.732050808
1.732050808
2
2
2.236067977
2.236067977
2.236067977
2.449489743
2.236067977
2.236067977
2.449489743
2.828427125
1.732050808
2.449489743
2.645751311
1.414213562
2.236067977
2
1.414213562
2.449489743
2.236067977
2.449489743
2.236067977
1.414213562
2.449489743
2.236067977
2.645751311
1.414213562
2.449489743
1.414213562
1.732050808
2
2
2.236067977
2.236067977
2.236067977
2.449489743
2.236067977
2
1.732050808
2.236067977
2.236067977
2
2.236067977
2
1
1.732050808
2
2
2.236067977
2.236067977
1.732050808
2
2.236067977
2.236067977
2
2.449489743
1.732050808
2.449489743
1.732050808
2
2.828427125
2.449489743
1.732050808
2.449489743
2.645751311
2.449489743
2.449489743
2.645751311
2.645751311
2.449489743
1.414213562
2.449489743
2.645751311
2.449489743
2.449489743
2.645751311
2.236067977
2.236067977
1.414213562
1.732050808
2.645751311
2.236067977
1.414213562
2.645751311
2
2.645751311
1.732050808
2.449489743
2
2.645751311
1
2.645751311
2.828427125
3
2.645751311
2.449489743
2
2.236067977
2.449489743
1.414213562
2.828427125
1.732050808
2.828427125
2.236067977
2
2.449489743
2.236067977
2.645751311
1.414213562
2.449489743
2
2.236067977
2.449489743
2.449489743
2.645751311
1.732050808
1
2.645751311
1.732050808
2
2.236067977
2
2.236067977
1.732050808
2.449489743
2
2.645751311
1
2.236067977
2.449489743
2.645751311
2.236067977
2
2.449489743
2.449489743
1.414213562
2.236067977
2
1.732050808
2
1.414213562
2.236067977
1.732050808
2.449489743
1.414213562
2
2.236067977
2.449489743
2
1.732050808
2.645751311
2.828427125
2.236067977
2.828427125
2.236067977
2
2.449489743
2.236067977
2.645751311
1.414213562
2.828427125
1.414213562
2.236067977
2.449489743
2.449489743
2.645751311
2.236067977
2.645751311
1.414213562
1.732050808
2.449489743
2
2.236067977
1.732050808
2.828427125
2
2.449489743
1.732050808
2
1.414213562
1
2.645751311
2.645751311
2
2.236067977
2.236067977
2.449489743
2.449489743
2.236067977
1.732050808
2.236067977
2.449489743
2.645751311
2.645751311
2.449489743
1.414213562
1.732050808
2.449489743
2.449489743
2.645751311
2.236067977
2.449489743
2.449489743
2.449489743
1.732050808
1.414213562
2
1
2.645751311
2.236067977
1.732050808
2
1.414213562
2.236067977
2.236067977
2.236067977
2
2.236067977
2.236067977
1.414213562
2.449489743
2
1.732050808
2.236067977
2
2.449489743
1.414213562
1.732050808
2
2
2.236067977
2.236067977
1.732050808
1
2.449489743
2
2.449489743
2.645751311
2.828427125
2.449489743
2.236067977
2.645751311
1.414213562
2.236067977
2.645751311
2.236067977
2
2.236067977
1.732050808
2.449489743
2.449489743
2.645751311
2.236067977
2.645751311
2.449489743
2.645751311
2.236067977
2
2.828427125
2.828427125
2
2.236067977
2
2.449489743
2.645751311
2.236067977
2.449489743
2.645751311
2.828427125
2.449489743
2.236067977
2.236067977
1.732050808
2
2
2.236067977
2.236067977
1
1
1.414213562
2
1.414213562
1.732050808
2.236067977
1.732050808
2.236067977
2.449489743];



% add in gibson shape differences
tempsimmeasures = zeros(325, size(simmeasures,2)+1);
tempsimmeasures(:,end)=simmeasures(:,end);
tempsimmeasures(:,1:end-1) = simmeasures;

tempsimmeasures(:,end-1) = gibsonshapemeasure;


simmeasures = tempsimmeasures;

% had to strip the column headers which are here
% last column is from Eagleman database

columnheaders = {'BC1989 SimRate Upper'	'BC1989 SimRate Lower'	'PG 1979 DiffRate Avg' ...
    'PG1979 RT Avg'	'KB FeatDist'	'g1979Confus'	'ggm1983ConfusA' ...
    'ggm1983ConfusB'	'Popp1964ConfusLower'	'BR1969ConfusLower' ...
    'BR1968ConfusUpper'	'Lewand 2000 FreqDiff'	'Lewand 2000 FreqSum' ...
    'Lewand FreqDiff/FreqSum'	'Ordinality Diff'	'OrdinalitySum'	'Ordinality Diff/Sum' ...
    'Gibson shape','EaglemanRGBdists'};


% can just plot them all

for i=1:size(simmeasures,2)-1
    
    
%     compute a correlation between the measures
%  linear regression
    p=polyfit(simmeasures(:,i),simmeasures(:,end),1);
%   points predicted by fit line
    yfit=polyval(p,simmeasures(:,i));
%   get the residuals
    yresid=simmeasures(:,end)-yfit;
%    sumsquaredresidual
    SSresid=sum(yresid.^2);
%     total sum of squares
    SSTotal=(length(simmeasures(:,end))-1)*var(simmeasures(:,end));
%     rsquared
    rsq = 1 - SSresid/SSTotal;
%     

    figure('Name',columnheaders{i},'Color',[1 1 1]);
         plot(simmeasures(:,i),yfit,'co','MarkerSize',10,'MarkerFaceColor','c');
         hold on;
    plot(simmeasures(:,i),simmeasures(:,end),'r.','MarkerSize',1);
    xlabel(columnheaders{i});
    ylabel('mean rgb distance');
    box off;
     text(simmeasures(:,i),simmeasures(:,end),letterpairs);
     hold on;
     title([num2str(rsq) '% variance explained ' 'r=' num2str(sqrt(rsq))]);

end




% correlation of measures with one another
figure;
plotmatrix(simmeasures)


[rho pval] = corr(simmeasures,'rows','pairwise','type','Spearman');

figure('Name','Correlation between lettersimilarity variables','Color',[1 1 1]);

subplot(1,2,1);
imagesc(rho);
colorbar;

subplot(1,2,2);
imagesc(pval);
colorbar;




