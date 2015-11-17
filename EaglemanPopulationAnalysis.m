
% wrapper for analyzing the eagleman database generally

% data and functions on my laptop


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
letters = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

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

% get LAB
p_Lab = RGB2Lab(p_rgb);






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
% also rounding is a bad idea.  if data are random uniform distributed when
% rounded they will appear gaussian as size of space which rounds to a
% point is reduced at edges.  so for example if rounding points are 0 .5 1,
% then 0-.25 rounds to 0, .25 to .75 rounds to .5, and .75 to 1 rounds to
% 1.  the plots as is actually underestimate the size of the effect because
% the rounding is working against it
% easiest is probably to divide rgb into cubes and then interpolate to the
% center of each cube.

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
% rgb_dists = nan(325,6588);
lab_dists = nan(325,6588);

% % calculate distance matrix for each subject
% for i=1:length(p_rgb)
% %     fill with calculated pdists
% % comparisons with nans appear as nans
%     rgb_dists(:,i)=pdist(squeeze(p_rgb(i,:,:)));
% 
% end

% calculate distance matrix for each subject
for i=1:length(p_Lab)
%     fill with calculated pdists
% comparisons with nans appear as nans
    lab_dists(:,i)=pdist(squeeze(p_Lab(i,:,:)));

end

% visualize average distance matrix

% meanrgbdist = nanmean(rgb_dists,2);
% try median
meanlabdist = nanmedian(lab_dists,2);

 figure('Name','histogram of average distance between pairs of letters in rgb','Color',[1 1 1]);
 
 hist(meanlabdist)
 xlabel('euclidian distance in L*a*b*');
 ylabel('count of letter pairs');
 
%  
% %  let's get histogram of distances for every possible letter pair
% 
% if ~exist('pairwisedists','dir')
%     mkdir('pairwisedists');
% end
% 
% for i=1:325
%     
%  
%  figure('Name',['histogram of L*a*b* distances for ' letterpairs{i}],'Color',[1 1 1]);
%  
%  hist(lab_dists(i,:));
%  xlabel('euclidian distance in L*a*b*');
%  ylabel('num of subjects');
%  box off;
%  
% %  save the figure
% 
% saveas(gcf,['pairwisedists/letterpairs{i} 'hist.png'],'png');
% 
% 
% % close
%     close(gcf);
%  
% end
%  
 

% what would the distances look like if they were chosen randomly?

% from RGB2Lab.m
%  Values for L are in the
% range [0,100] while a and b are roughly in the range [-110,110].  The
% output is of type double.

% generate a dataset

% rnadom values on unit interval
randLab = rand(size(p_Lab));
% scale to lab
% L*
randLab(:,:,1) = 0+(100-0).*randLab(:,:,1);
% a*
randLab(:,:,2) = (110--110).*randLab(:,:,2)-110;
% b*
randLab(:,:,3) = (110--110).*randLab(:,:,3)-110;

% randLabdists
% calculate distance matrix for each subject
randlab_dists = nan(325,6588);

for i=1:length(randLab)
%     fill with calculated pdists
% comparisons with nans appear as nans
    randlab_dists(:,i)=pdist(squeeze(randLab(i,:,:)));

end

% try median
medrandlabdist = nanmedian(randlab_dists,2);

 figure('Name','histogram of average distance between pairs of letters in L*a*b* random data','Color',[1 1 1]);
 
 hist(medrandlabdist)
 xlabel('euclidian distance in L*a*b*');
 ylabel('count of letter pairs');
 box off;

 
 
%  random data looks very different.  in fact average distance would be
%  much larger and much more tightly distributed.  however, this is
%  probably because the gamut of Lab is larger than RGB so its not
%  surprising that the distances are smaller
% let's make random rgb data, pass it through RGB2Lab and see what we get

 
 
% random values on unit interval
randrgb = rand(size(p_rgb));
% convert to lab
randrgb2lab = RGB2Lab(randrgb);

% its obvious that RGB and LAB have different ranges.  for example the
% minimum L value after suhuffling is around 10 instead of 0.  the max is
% >99
% when the data are randomly chosen in RGB first it now appears as if the
% average distance between a pair of letter is larger than if you just
% chose points randomly.  this is probably not that surprising since there
% is a tendencry to pick extreme values in at least one of the RGB
% channels.


% randrgbdists
% calculate distance matrix for each subject
randrgb2lab_dists = nan(325,6588);

for i=1:length(randrgb2lab)
%     fill with calculated pdists
% comparisons with nans appear as nans
    randrgb2lab_dists(:,i)=pdist(squeeze(randrgb2lab(i,:,:)));

end

% try median
medrandrgb2labdist = nanmedian(randrgb2lab_dists,2);

 figure('Name','histogram of average distance between pairs of letters in L*a*b* random data in RGB transformed','Color',[1 1 1]);
 
 hist(medrandrgb2labdist)
 xlabel('euclidian distance in L*a*b*');
 ylabel('count of letter pairs');
 box off;
 
 
%  plot against linear distance in sequence

% 
% this will be handy for other things
letterseqdists = pdist((1:26)');


% figure('Name','L*a*b* distance vs alphabet distance','Color',[1 1 1]);
% 
% plot(letterseqdists,meanlabdist,'ro');
% box off;
% xlabel('alphabetic distance');
% ylabel('mean lab distance');

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
% make a new matrix and add one extra column
tempsimmeasures = zeros(325, size(simmeasures,2)+1);

tempsimmeasures(:,1:end-1) = simmeasures;

tempsimmeasures(:,end-1) = gibsonshapemeasure;


simmeasures = tempsimmeasures;

% make the last column the lab measures instead
simmeasures(:,end) = meanlabdist;

% had to strip the column headers which are here
% last column is from Eagleman database

columnheaders = {'BC1989 SimRate Upper'	'BC1989 SimRate Lower'	'PG 1979 DiffRate Avg' ...
    'PG1979 RT Avg'	'KB FeatDist'	'g1979Confus'	'ggm1983ConfusA' ...
    'ggm1983ConfusB'	'Popp1964ConfusLower'	'BR1969ConfusLower' ...
    'BR1968ConfusUpper'	'Lewand 2000 FreqDiff'	'Lewand 2000 FreqSum' ...
    'Lewand FreqDiff div FreqSum'	'Ordinality Diff'	'OrdinalitySum'	'Ordinality Diff div Sum' ...
    'Gibson shape','EaglemanLabdists'};


% can just plot them all


if ~exist('correlationsAveAcrossSubs','dir')
    mkdir('correlationsAveAcrossSubs');
end

for i=1:(size(simmeasures,2)-1)
    
    
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
%   compare to matlab output
    [s_rho, pval] = corr(simmeasures(:,i),simmeasures(:,end),'type','Spearman');

    figure('Name',columnheaders{i},'Color',[1 1 1]);
         plot(simmeasures(:,i),yfit,'co','MarkerSize',10,'MarkerFaceColor','c');
         hold on;
    plot(simmeasures(:,i),simmeasures(:,end),'r.','MarkerSize',1);
    xlabel(columnheaders{i});
    ylabel('mean L*a*b* distance');
    box off;
     text(simmeasures(:,i),simmeasures(:,end),letterpairs);
     hold on;
     title([num2str(rsq) '% variance explained ' 'r=' num2str(sqrt(rsq)) ' : '...
         'rho = ' num2str(s_rho) ' p = ' num2str(pval)]);
     
     saveas(gcf,['correlationsAveAcrossSubs/' columnheaders{i} '.png'],'png');
     close gcf;

end





% 
% 
% % correlation of measures with one another
% figure;
% plotmatrix(simmeasures)


[rho pval] = corr(simmeasures,'rows','pairwise','type','Spearman');

figure('Name','Correlation between lettersimilarity variables','Color',[1 1 1]);

subplot(1,2,1);
imagesc(rho);
colorbar;

subplot(1,2,2);
imagesc(pval);
colorbar;



% in both the brang and watson papers they bin the letter pairs into groups
% based on the average similarity in color distance
% simmeasures is the average distance in RGB (in future we will match the
% color space as well)

[sortedsim, sortedindx] = sort(simmeasures(:,end));

% average our responses

% since the inflation of the correlation is dependent on the smoothing we
% can play with this
binsize = 5;

% reshape for averaging
sortedsim=reshape(sortedsim,binsize,325/binsize);
% average
ave_sortedsim=mean(sortedsim);

% sort our letter pairs to match
sortedletterpairs = letterpairs(sortedindx);

% if ~exist('correlationsAveAcrossSubsAndLetterPairs','dir')
%     mkdir('correlationsAveAcrossSubsAndLetterPairs');
% end
% 
% 
% % fastest way is to reshape each matrix and compute average?
% 
% for i=1:size(simmeasures,2)-1
%     
%     % ave sorted predictor
%     s=simmeasures(sortedindx,i);
%     % reshape for averaging (checked this bookkeeping and it looks right)
%     s=reshape(s,binsize,325/binsize);
%     % average the columns
%     ave_s=mean(s);
%     
%     %     compute a correlation between the measures
%     %  linear regression
%     p=polyfit(ave_s,ave_sortedsim,1);       
%     %   points predicted by fit line
%     yfit=polyval(p,ave_s);
%     %   get the residuals
%     yresid=ave_sortedsim-yfit;
%     %    sumsquaredresidual
%     SSresid=sum(yresid.^2);
%     %     total sum of squares
%     SSTotal=(length(ave_sortedsim)-1)*var(ave_sortedsim);
%     %     rsquared
%     rsq = 1 - SSresid/SSTotal;
%     %
%    %   compare to matlab output
%     [s_rho, pval] = corr(ave_s',ave_sortedsim','type','Spearman');
% 
%     
%     figure('Name',columnheaders{i},'Color',[1 1 1]);
%     plot(ave_s,yfit,'co','MarkerSize',10,'MarkerFaceColor','c');
%     hold on;
%     plot(ave_s,ave_sortedsim,'r.','MarkerSize',10);
%     xlabel(columnheaders{i});
%     ylabel('mean L*a*b* distance');
%     box off;
% %     text(simmeasures(:,i),simmeasures(:,end),letterpairs);
%     hold on;
%       title([num2str(rsq) '% variance explained ' 'r=' num2str(sqrt(rsq)) ' : '...
%          'rho = ' num2str(s_rho) ' p = ' num2str(pval)]);
%      
%        saveas(gcf,['correlationsAveAcrossSubsAndLetterPairs/' columnheaders{i} '.png'],'png');
%      close gcf;
% 
% end








% so this is really the wrong approach.  instead we should build a full
% regression model with subject effects.  the one difficulty is that about
% 35% of the subjects failed to match 1 or more letters.  however this
% still leaves ~4300 subjects which is not bad.  we could also try the lme
% approach but I'm not sure how much it would help us given the level of
% complexity it introduces.

% there ia another layer of difficulty in assessing the correlations.   the
% assumption is that all the samples are independent, but since they are
% all pairwise distances, that is not really true. for exaple, the measures
% AB, AC, ... AZ all measure distances from A.  if A is far from everything
% else, than these will tend to be large.  Waskom said people are aware of
% this problem when doing RSM analyses generally (namely that your null
% hypothesis might not really be 0).  They use the following method,
% compute the correlation for each subject.  The correlation can be bootstrapped across
% subjects (so random fx if you need).

% for each subject you can also compute a null distribution using a
% permutation approach.  this is done by shuffling the rows and columns of
% the matrix while retaining their connection to one another.  so for
% example  rowA columnA   could go to rowB columnB.  this mixes up the
% distances while maintaining the dependences that arise from having one
% endpoint in the distance across 25 comparisons....what happens with the
% diagonal?



% can plot distributions correlations of individual subjects with each measures

for i=1:(size(simmeasures,2)-1)
    
    

%   compare to matlab output
    [rho, pval] = corr(simmeasures(:,i),lab_dists,'type','Spearman','row','pairwise');

    figure('Name',columnheaders{i},'Color',[1 1 1]);
    subplot(1,2,1);
       hist(rho,-1:.1:1);
    xlabel('rho');
    ylabel('count of subjects');
    box off;
     hold on;
     set(gca,'XLim',[-1 1]);
     subplot(1,2,2);
       hist(pval,0:.025:1);
    xlabel('pval');
    ylabel('count of subjects');
    box off;
     hold on;
     set(gca,'XLim',[0 1]);
  

end




% % % % % % % % % % % % % % % % % % % % 

% % %   made these figures, might have to rerun to do any stats

% % % % % % % % % % % % % % % % % % % % % % % % % % % % 




% this still isn't right as the correlations are still calculated as though
% all the points are independent and compared to a null which is
% distributed around 0. one way to approach this is to summarize each
% subject against its own permuted null.  

% 1.  for each subject we have a rho
% 2.  for each subject permute data and calculate rho many times
% 3.  percentile of data rho within permuted rho distribution


% % distancemeasure = 13; %pick a column to do
% 
% numpermutes = 100;
% 
% permute_rhos = zeros(numpermutes,1); %for each subjects null
% 
% pctlarger = nan(length(lab_dists),1); %percent of shuffles larger than observed for each subject
% 
% if ~exist('bootstrappedInd','dir')
%     mkdir('bootstrappedInd');
% end
% 
% for distancemeasure=1:(length(columnheaders)-1)
%     distancemeasure
% % for each subject
% for s=1:length(lab_dists)
% %     s
% %     
%     
% %     there must be a faster way
% % 1. make distance vector for a subject into full matrix
%     labsquare = squareform(lab_dists(:,s));
%     
% %     do permutations
%     for p= 1:numpermutes
% % 2. get permutationindex
%     [pletters pindex] = shuffle(1:26);
% % 3. permute cols
%     plabsquare = labsquare(:,pletters);
% % 4. permute rows
%     plabsquare =plabsquare(pletters,:);
% % 5. unpack upper diagonal
% % this cheat won't work if some distances are actually zero
%     tplabsquare = plabsquare';
%     permdist = tplabsquare(tril(true(size(plabsquare)),-1));
%     
%     
% % 6.compute a correlation
%     permute_rhos(p) = corr(simmeasures(:,distancemeasure),permdist, ...
%         'type','Spearman','rows','pairwise');
%     end    
% 
% % 7.  get actual correlation
%     rho = corr(simmeasures(:,distancemeasure),lab_dists(:,s), ...
%         'type','Spearman','rows','pairwise');
% % 8. get percent of correlations in permute larger than observed
%     pctlarger(s) = length(find(permute_rhos>=rho))/numpermutes;
%     
%     
%     
% end
% 
% 
% 
% % make a figure
% figure('Name',columnheaders{distancemeasure} , 'Color', [1 1 1])
% hist(pctlarger,0:.05:1);
% box off;
% ylabel('count of subjects');
% xlabel('percentile of observed subject correlation in distribution of nulls');
% 
% saveas(gcf,['bootstrappedInd/' columnheaders{distancemeasure}, '.png'],'png');
% close gcf;
% end
% 
