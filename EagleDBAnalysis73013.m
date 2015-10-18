

% it was then updated with some user ids

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

% some usefuls stuff
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




% show me the nans in the database
s_nansdistribution


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





% visualize the whole database
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

%% plot colors for matches that are close to magnet set
% % % whole data set with modes for each letter


figure('name', 'whole database', 'Color', [1 1 1],'Position',get(0,'ScreenSize'));

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:length(syntype)/10,:,:);
theResult = rgb.eagle;

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
box off;


% let's plot each group

% magnets
% matches
figure('name', [ num2str(length(find(syntype==2))) ...
    ' subjects with more than ' num2str(magnetthreshold) ...
    ' matches to letter set'], 'Color', [1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,3,1);

% make a graphical legend 10% as tall as the result
theLegend = rgb.magnets(1:length(find(syntype==2))/10,:,:);
theResult = rgb.eagle(find(syntype==2),:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title([ num2str(length(find(syntype==2))) ...
    ' subjects with more than ' num2str(magnetthreshold) ...
    ' matches to letter set'])

% if you wanted to compare it to the earlier figure, you would shuffle the
% data instead of the magnet set and then sort.  the magnet set is just one
% throw instead of n?
% would be better to use the rich, but don't have rgb values, just
% simulated labels

% All the shuffled matches
% figure('name', ['top ' num2str(length(find(syntype==2))) ' shuffled matches to letter set'], 'Color', [1 1 1]);

subplot(1,3,3);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% this is now
[Y, ranking] = sort(nummatches.eagleShuffleByCol, 'descend');

best = ranking(1:length(find(syntype==2)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.magnets(1:round(length(find(syntype==2))/10),:,:);
theResult = rgb.eagleShuffledByCol(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title('All col shuffled matches to magnet set')


subplot(1,3,2);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% this is now
[Y, ranking] = sort(nummatches.eagleShuffleByRow, 'descend');

best = ranking(1:length(find(syntype==2)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.magnets(1:round(length(find(syntype==2))/10),:,:);
theResult = rgb.eagleShuffledByRow(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title('All row shuffled matches to magnet set')
% 
% % makes this into a nice sfn figure
% fsize=get(gcf,'Position');
% set(gcf,'Position',[5 5 5*fsize(3) 5*fsize(4)]);
% saveas(gcf,'magnetsAndShuffledRGB.jpg','jpg');
% 



% how are nans distributed across letters for this subset?
magsubindx = find(syntype==2);
figure('Name','distribution of nans across letters','Color',[1 1 1]);

subplot(1,2,1);
% figure where nans are black and color matches are white
% get one rgb column of data for subs x letters

imagesc(nansmatrix(magsubindx,:));  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(nansmatrix(magsubindx,:))/length(magsubindx));;
box off;
xlabel('letters');
ylabel('number of times not matched');
set(gca,'XTick',[1:26],'XTickLabel',letters);







% histogram
% colors come from fp_visualizeEagDB
makeLetterXColorHist(dbNumbered(find(syntype==2),:));
set(gcf,'Name','just magnet syns');





% want to slip a table in here too in the command window
disp('magnet proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==2),i),11);
    counts  = counts/length(find(syntype==2));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end


% high fq
figure('name', [ num2str(length(find(syntype==1))) ...
    ' subjects with more than ' num2str(fqthreshold) ...
    ' matches to high fq template'], 'Color', [1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,3,1);

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:length(find(syntype==1))/10,:,:);
theResult = rgb.eagle(find(syntype==1),:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title([ num2str(length(find(syntype==1))) ...
    ' subjects with more than ' num2str(magnetthreshold) ...
    ' matches to high fq template']);




% if you wanted to compare it to the earlier figure, you would shuffle the
% data instead of the magnet set and then sort.  the magnet set is just one
% throw instead of n?  would be good to use rich, but it is only simulated
% labels, not really good to simulate rgb values, though we could

% All the shuffled matches
% figure('name', ['top ' num2str(length(find(syntype==1))) ' shuffled matches to high fq template'], 'Color', [1 1 1]);
subplot(1,3,3);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% is now
[Y, ranking] = sort(nummatches.eagleShuffleByCol, 'descend');
best = ranking(1:length(find(syntype==1)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:round(length(find(syntype==1))/10),:,:);
theResult = rgb.eagleShuffledByCol(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18);
title(['top ' num2str(length(find(syntype==1))) ' col shuffled matches to high fq template']);

subplot(1,3,2);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% is now
[Y, ranking] = sort(nummatches.eagleShuffleByRow, 'descend');
best = ranking(1:length(find(syntype==1)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:round(length(find(syntype==1))/10),:,:);
theResult = rgb.eagleShuffledByRow(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18);
title(['top ' num2str(length(find(syntype==1))) ' row shuffled matches to high fq template']);

% 
% % makes this into a nice sfn figure
% fsize=get(gcf,'Position');
% set(gcf,'Position',[5 5 5*fsize(3) 5*fsize(4)]);
% saveas(gcf,'svgplots/SFN2013AndShuffledRGB.jpg','jpg');
% 



% how are nans distributed across letters for this subset?
fqsubindx = find(syntype==1);
figure('Name','distribution of nans across letters for fq syns','Color',[1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,2,1);
% figure where nans are black and color matches are white
% get one rgb column of data for subs x letters

imagesc(nansmatrix(fqsubindx,:));  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(nansmatrix(fqsubindx,:))/length(fqsubindx));
box off;
xlabel('letters');
ylabel('number of times not matched');
set(gca,'XTick',1:26,'XTickLabel',letters);










% histogram
% colors come from fp_visualizeEagDB
makeLetterXColorHist(dbNumbered(find(syntype==1),:));
set(gcf,'Name','just high fq syns');

% want to slip a table in here too in the command window
disp('high fq proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==1),i),11);
    counts  = counts/length(find(syntype==1));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % NO GROUP MATCHES AND SHUFFLING

% remainder
figure('name', [ num2str(length(find(syntype==0))) ...
    ' subjects with no group'], 'Color', [1 1 1]);

subplot(1,3,1);
% make a graphical legend 10% as tall as the result
% theLegend = rgb.fq(1:length(find(syntype==2))/10,:,:);
theResult = rgb.eagle(find(syntype==0),:,:);
theLegend = rgb.fq(1:round(length(find(syntype==0))/10),:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title([ num2str(length(find(syntype==0))) ...
    ' subjects with no group']);

% All the shuffled matches
% figure('name', ['top ' num2str(length(find(syntype==1))) ' shuffled matches to high fq template'], 'Color', [1 1 1]);
subplot(1,3,3);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% is now
[Y, ranking] = sort(nummatches.eagleShuffleByCol, 'descend');
best = ranking(1:length(find(syntype==0)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:round(length(find(syntype==0))/10),:,:);
theResult = rgb.eagleShuffledByCol(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18);
title(['top ' num2str(length(find(syntype==1))) ' col shuffled matches to high fq template']);

subplot(1,3,2);

% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% is now
[Y, ranking] = sort(nummatches.eagleShuffleByRow, 'descend');
best = ranking(1:length(find(syntype==0)));

% make a graphical legend 10% as tall as the result
theLegend = rgb.fq(1:round(length(find(syntype==0))/10),:,:);
theResult = rgb.eagleShuffledByRow(best,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18);
title(['top ' num2str(length(find(syntype==1))) ' row shuffled matches to high fq template']);






% how are nans distributed across letters for this subset?
ngsubindx = find(syntype==0);
figure('Name','distribution of nans across letters for no group syns','Color',[1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,2,1);
% figure where nans are black and color matches are white
% get one rgb column of data for subs x letters

imagesc(nansmatrix(ngsubindx,:));  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(nansmatrix(ngsubindx,:))/length(ngsubindx));;
box off;
xlabel('letters');
ylabel('number of times not matched');
set(gca,'XTick',[1:26],'XTickLabel',letters);





% histogram
% colors come from fp_visualizeEagDB

% histogram
% colors come from fp_visualizeEagDB
makeLetterXColorHist(dbNumbered(find(syntype==0),:));
set(gcf,'Name','no group syns');



% do a histogram of just the labels that don't come from the magnet
% synesthetes
makeLetterXColorHist(dbNumbered(find(syntype~=1),:));
set(gcf,'Name','non magnet syns');



format short g;

% want to slip a table in here too in the command window
disp('no group proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==0),i),11);
    counts  = counts/length(find(syntype==0));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end


% want to slip a table in here too in the command window
disp('magnet proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==2),i),11);
    counts  = counts/length(find(syntype==2));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end

% want to slip a table in here too in the command window
disp('high fq proportions');
disp(['letter  black  white  red    green  yellow blue   brown  purple pink   orange grey']);
for i=1:26
    [counts, bins]=hist(dbNumbered(find(syntype==1),i),11);
    counts  = counts/length(find(syntype==1));
    disp([letters(i) '       ' num2str(counts,'%0.2f   ')]);
    
end









% correlation of nans among the subgroups

figure('Name','Correlation of non matches to letters across groups','Color',[1 1 1],'Position',get(0,'ScreenSize'));

% values to correlate
fqnans=sum(nansmatrix(fqsubindx,:)/length(fqsubindx));
magnans=sum(nansmatrix(magsubindx,:)/length(magsubindx));
ngnans=sum(nansmatrix(ngsubindx,:)/length(ngsubindx));

subplot(1,3,1);
scatter(fqnans,magnans);
text(fqnans,magnans,letters,'FontSize',14);
hold on;
plot(0:.01:.2,0:.01:.2,'k--');
set(gca,'YLim',[0 .2],'XLim',[0 .2],'YTick',0:.05:.2,'XTick',0:.05:.2);
xlabel('proportion nans in fq group');
ylabel('proportion nans in magnet group');
box off;
axis square;


subplot(1,3,2);
scatter(fqnans,ngnans);
text(fqnans,ngnans,letters,'FontSize',14);
hold on;
plot(0:.01:.2,0:.01:.2,'k--');
set(gca,'YLim',[0 .2],'XLim',[0 .2],'YTick',0:.05:.2,'XTick',0:.05:.2);
xlabel('proportion nans in fq group');
ylabel('proportion nans in non group');
box off;
axis square;

subplot(1,3,3);
scatter(ngnans,magnans);
hold on;
plot(0:.01:.2,0:.01:.2,'k--');
text(ngnans,magnans,letters,'FontSize',14);
set(gca,'YLim',[0 .2],'XLim',[0 .2],'YTick',0:.05:.2,'XTick',0:.05:.2);
xlabel('proportion nans in no group');
ylabel('proportion nans in magnet group');
box off;axis square;



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % %   do age analysis

load('AgeMatrix.mat');
% loads variable called AgeMatrix
% columns are subjectid batteryid YOB agein2014

% should be able to search through this using USERID which i think
% corresponds to batterid for the subset of our subjects with age data
% userid is the same length as the data so has the same index

for i=1:length(userid)
    %     see if gc synesthete gave a usable age
    indx = find(AgeMatrix(:,2)==userid(i));
    %     if not then record a nan for age and dob
    if isempty(indx)
        subage(i) = nan;
        subdob(i) = nan;
    else
        %         record 2014 age and dob
        subage(i) = AgeMatrix(indx,4);
        subdob(i) = AgeMatrix(indx,3);
    end
end


% make some figures
% histograms of subjects by year, fq subs by year, and mag subs by year
figure('name','subject date of birth data','Color', [1 1 1],'Position',get(0,'ScreenSize'))
% allsubjects
subplot(1,3,1);
hist(subdob,1930:5:2010);
box off;
ylabel('num subjects');
xlabel('year born');
title('all subjects');
set(gca,'YLim',[0 1800]);

% fq subjects
subplot(1,3,2);
hist(subdob(syntype==1),1930:5:2010)
box off;
ylabel('num subjects');
xlabel('year born');
title('fq subjects');
set(gca,'YLim',[0 1800]);

% magnet subjects
subplot(1,3,3);
hist(subdob(syntype==2),1930:5:2010)
box off;
ylabel('num subjects');
xlabel('year born');
title('magnet subjects');
set(gca,'YLim',[0 1800]);

% analyse the proportion of the database having magnets overtime

magdob = subdob(syntype==2);
notmagdob = subdob(syntype~=2);
fqdob = subdob(syntype==1);

% how many are there?
nummag = sum(~isnan(magdob));
numnotmag = sum(~isnan(notmagdob));


% histogram data in 5 year intervals
bins = 1925:5:2005;

maghist = histc(magdob,bins);
notmaghist = histc(notmagdob,bins);
allhist = histc(subdob,bins);
fqhist = histc(fqdob,bins);

figure('name','prevalance of magnet synesthetes over time','Color',[1 1 1],'Position',get(0,'ScreenSize'));
plot(bins,maghist./allhist,'ro-','MarkerFaceColor',[1 0 0]);
box off;
xlabel('subjects born in 5 years after...');
ylabel('proportion of subjects >=10 matches to magnets');

plot2svg('prevalenceOverTimeFull.svg',gcf,'svg');


% this might look better as a bar chart
figure('name','prevalance of magnet synesthetes over time','Color',[1 1 1],'Position',get(0,'ScreenSize'));
bar(bins,maghist./allhist,1);
box off;
set(gca,'XTick',bins);
xlabel('year subjects with 10 or more matches were born');
ylabel('proportion of subjects >=10 matches to magnets');

plot2svg('prevalenceOverTimeFullBars.svg',gcf,'svg');





figure('name','prevalance of fq and mag synesthetes over time','Color',[1 1 1],'Position',get(0,'ScreenSize'));
plot(bins,maghist./allhist,'ro-','MarkerFaceColor',[1 0 0]);
hold on;
plot(bins,fqhist./allhist,'ko-','MarkerFaceColor',[0 0 0]);
box off;
xlabel('subjects born in 5 years after...');
ylabel('proportion of subjects >=10 matches to magnets');

plot2svg('prevalenceOverFQvsMAG.svg',gcf,'svg');



% suppose we want confidence intervals on the histograms.  its weird to say
% that because we only have the one measure, but very different amounts of
% data to work with over the years (very few subjects in their 80s or under
% the age of 10) many in the 30-50 range

% [bins' allhist' fqhist' maghist']

% 
%         1925           0           0           0
%         1930           6           0           0
%         1935           8           0           0
%         1940          27           3           0
%         1945          42           6           0
%         1950          71          10           0
%         1955         104          19           0
%         1960         126           8           0
%         1965         186          27          13
%         1970         293          24          41
%         1975         535          76          80
%         1980         981         142         126
%         1985        1460         284          77
%         1990        1632         389          36
%         1995         607         147           8
%         2000         128          23           4
%         2005          14           2           1

% so the algorithm is to randomly sample with replacement from our
% distribution many times and then generate statistics across the bins
% so our subject ages are in subdob
    nbstraps = 1000;

    fqstraps = [];
    mgstraps = [];

    % histogram data in 5 year intervals
    bins = 1925:5:2005;

    % probably doesn't need to be a for loop?
    % for each bootstrap
    for i=1:nbstraps
        % get index to our random sample
        rsample =  randi(length(subdob),[length(subdob),1]);
    %    get type of synesthete for sample index
        rsampsyntype = syntype(rsample);
    %     get dobs for sample index
        rsampsubdob = subdob(rsample)';
        %     turn these into binned data for different types of synesthetes
        rsampallhist = histc(rsampsubdob,bins);
        rsampfqhist = histc(rsampsubdob(rsampsyntype==1),bins);
        rsampmaghist = histc(rsampsubdob(rsampsyntype==2),bins);
    %     now into proportions which we will store
        fqstraps(i,:) = rsampfqhist./rsampallhist;
        mgstraps(i,:) = rsampmaghist./rsampallhist;


    end

% figure with all our bootstraps
    figure('Name','All bootstraps of fq and mag','Color',[1 1 1],'Position',get(0,'ScreenSize'));
    plot(bins(1:15), fqstraps(:,1:15),'k--');
    hold on;
    plot(bins(1:15), mgstraps(:,1:15), 'r--');
    xlabel('year born');
    ylabel('perecent of population');


% figure with shaded 95% confidence intervals
    
%  get intervals

fqci = prctile(fqstraps,[2.5 97.5]);
mgci = prctile(mgstraps,[2.5 97.5]);

% errorbars use differences from mean not actual values
fqmed = nanmedian(fqstraps);%median bootstrapped high frequency syns
fqer(1,:) = abs(fqmed - fqci(1,:));
fqer(2,:) = abs(fqmed + fqci(2,:));
mgmed = nanmedian(mgstraps);%median bootstrapped magnet syns
mger(1,:) = abs(mgmed - mgci(1,:));
mger(2,:) = abs(mgmed + mgci(2,:));




figure('Name','95% ci of bootstraps of fq and mag','Color',[1 1 1],'Position',get(0,'ScreenSize'));

whichbins = [3:15];
errorbar3(bins(whichbins),fqmed(whichbins),fqci(:,whichbins),1,'k');
hold on;
plot(bins(whichbins),fqmed(whichbins),'k');
errorbar3(bins(whichbins),mgmed(whichbins),mgci(:,whichbins),1,'r');
plot(bins(whichbins),mgmed(whichbins),'r');
axis on;
box on;
% grid on;

% another way to do this is to bootstrap the mean or confidence interval within
% each bin


% plot2svg('../PrevalenceOfLearning/JustTheMagnets/psychsci/prevalenceOverFQvsMAGbootstrapped.svg',gcf,'svg');





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % analysis of scores for different groups

LoadColorSequenceScores
% loads a matrix
%<pre>
% This export generated on 2014-03-06
% Fields are: users_id, BatteryId, LCAvgScore, NCAvgScore, GCAvgScore,
% WCAvgScore, MCAvgScore, SCNumCorrect, SCNumIncorrect, SCNumTotal, SCAccuracyPercent, SCMeanRTCorrect, SCMeanRTIncorrect, SCMeanRTTotal
% matrix has a lot of subjects, so once again we want to find ours. this time all of them should be in here somewhere       
% 50868          14
% oddly, batteryID again appears to correspond to the variable userid we
% already have


% matrix to hold our behavior
behavior = zeros(6588,14);

for i=1:length(userid)
%     find each subject
    indx = find(ColorSequenceScores(:,2)==userid(i));
%     add their behavior to our matrix
    behavior(i,:) = ColorSequenceScores(indx,:);

end



% color matching distance threshold is 1
% median dist for magnet syns
magcdist = median(behavior(find(syntype==2),3));
% 95% confidence interval for magnet syns
% [magcdistci, magcdiststat] = bootci(5000,{@nanmedian,behavior(find(syntype==2),3)});
magcdistci = bootci(5000,{@nanmedian,behavior(find(syntype==2),3)});
% everyone else
notmagcdist = median(behavior(find(syntype~=2),3));
notmagcdistci = bootci(5000,@median,behavior(find(syntype~=2),3));
% everyone
allcdist = median(behavior(:,3));
allcdistci = bootci(5000,@median,behavior(:,3));


% make plot and do a ttest
[hcdist pcdist pci pstats] = ttest2(behavior(find(syntype==2),3),behavior(find(syntype~=2),3))
figure('Name','magnet vs nonmagnet color matching scores','Color',[1 1 1],'Position',get(0,'ScreenSize'));
bar(1:2,[magcdist, notmagcdist]);
hold on;
box off;
% confidence intervals
errorbar2(1:2,[magcdist,notmagcdist],[magcdistci notmagcdistci],1,'k');
set(gca,'XTickLabel',{'magnet syns','rest of pop'});
ylabel('color matching score');
title(['t = ' num2str(pstats.tstat) ', p = ' num2str(pcdist)]);

% speeded classification accuracy
% median acc for magnets accuracy is column ll
magacc=median(behavior(find(syntype==2),11));
magaccci = bootci(5000,@median,behavior(find(syntype==2),11));
% median acc for not magnets
notmagacc=median(behavior(find(syntype~=2),11));
notmagaccci = bootci(5000,@median,behavior(find(syntype~=2),11));
% median acc for everything
allacc = median(behavior(:,11));
allaccci = bootci(5000,@median,behavior(:,11));


% make plot and do a ttest
[hcacc pcacc pci pstats] = ttest2(behavior(find(syntype==2),11),behavior(find(syntype~=2),11))
figure('Name','magnet vs nonmagnet color matching scores','Color',[1 1 1],'Position',get(0,'ScreenSize'));
bar(1:2,[magacc, notmagacc]);
hold on;
box off;
% confidence intervals
errorbar2(1:2,[magacc,notmagacc],[magaccci notmagaccci],1,'k');
set(gca,'XTickLabel',{'magnet syns','rest of pop'});
ylabel('speeded classification accuracy');
title(['t = ' num2str(pstats.tstat) ', p = ' num2str(pcacc)]);


% speeded classification correct rts
% median magnets
magrt=median(behavior(find(syntype==2),12));
% bootstrap 95% ci
magrtci = bootci(5000,@median,behavior(find(syntype==2),12));
% median not magnets
notmagrt=median(behavior(find(syntype~=2),12));
notmagrtci = bootci(5000,@median,behavior(find(syntype~=2),12));
% all subjects
allcrt=median(behavior(:,12));
allcrtci = bootci(5000,@median,behavior(:,12));


% make plot and do a ttest
[hcrt pcrt pci pstats] = ttest2(behavior(find(syntype==2),12),behavior(find(syntype~=2),12))
figure('Name','magnet vs nonmagnet color matching scores','Color',[1 1 1],'Position',get(0,'ScreenSize'));
bar(1:2,[magrt, notmagrt]);
hold on;
box off;
% confidence intervals
errorbar2(1:2,[magrt,notmagrt],[magrtci notmagrtci],1,'k');
set(gca,'XTickLabel',{'magnet syns','rest of pop'});
ylabel('speeded classification rt');
title(['t = ' num2str(pstats.tstat) ', p = ' num2str(pcrt)]);

% via kendrick's class
% - Randomization (or permutation) tests. Let's pose the null hypothesis
% that the two sets of data come from the same probability distribution
% (not necessarily Gaussian). Under the null hypothesis, the two sets of data
% are interchangeable, so if we aggregate the data points and randomly divide 
% the data points into two sets, then the results should be comparable to the 
% results obtained with the original data. So, the strategy is to generate
% random datasets, compute some statistic from these datasets (such as
% difference in means or difference in medians), and then compare the 
% resulting values to the statistic computed from the original data. We count
% the number of randomly obtained values that are more extreme than the actual
% observed value and divide this by the total number of simulations that were 
% run. The result is the p-value. Notice that we have used raw computational 
% power to calculate the p-value directly instead of relying on analytical 
% formulas (which are valid only if certain assumptions are met).

%
% so here we would randomly assign 400 subjects to 1 group and the
% remainder to another and find the difference in the median scores.  we
% create a distribution of these simulated medians and then find out where
% our actual data is on this distribution.  

% our observed difference in medians for color matching, acc, and rt
diff_cdistmedian = magcdist - notmagcdist;
diff_accmedian = magacc - notmagacc;
diff_rtmedian = magrt - notmagrt;

% now make a distribution of bootstrapped medians
nboots = 5000;
% variables holding bootstrapped distribution of medians
bs_diffcdistmedian = [];
bs_diffaccmedian = [];
bs_diffrtmedian = [];

% should be done as a matrix!!
for i=1:nboots
%     create a new sample by shuffling assignment of syntype
    shuffledsyntype = shuffle(syntype);
    
%     calculate new medians for cdist
    sh_magcdist = median(behavior(shuffledsyntype==2,3));
    sh_notmagcdist = median(behavior(shuffledsyntype~=2,3));
    bs_diffcdistmedian = [bs_diffcdistmedian sh_magcdist-sh_notmagcdist];
    
%     calculate new medians for accuracy
    sh_magacc = median(behavior(shuffledsyntype==2,11));
    sh_notmagacc = median(behavior(shuffledsyntype~=2,11));
    bs_diffaccmedian = [bs_diffaccmedian sh_magacc-sh_notmagacc];
    
    %     calculate new medians for rt
    sh_magrt = median(behavior(shuffledsyntype==2,12));
    sh_notmagrt = median(behavior(shuffledsyntype~=2,12));
    bs_diffrtmedian = [bs_diffrtmedian sh_magrt-sh_notmagrt];
    
end

figure('Name','bootstrapped medians for www.synesthete.org','Color',[1 1 1],'Position',get(0,'ScreenSize'));
% color matching
subplot(3,2,1);
bar(1:2,[magcdist, notmagcdist]);
hold on;
box off;
% confidence intervals
errorbar2(1:2,[magcdist,notmagcdist],[magcdistci notmagcdistci],1,'r');
set(gca,'XTickLabel',{'magnet syns','rest of pop'});
ylabel('rts for classification task');
% title();
subplot(3,2,2);
% calculate percentile less than observed value
p_rtmedian = 1-length(find(bs_diffrtmedian>diff_rtmedian))/nboots
% add mark for observed median
[counts bincenters] = hist(bs_diffrtmedian,30);
% turn it into a probability distribution
bar(bincenters,counts/nboots,1);
hold on;
% yvalue of correct bin from histogram
yval=counts(find(min(bincenters-diff_rtmedian)))/nboots;
% plot a line around that height of bin (will be small!)
plot(diff_rtmedian,yval-(yval/2):yval+(yval/2),'r','Linewidth',8);
box off;
xlabel('bootstrapped difference in median rts for classification');
ylabel('probability of difference in medians'); 
title(['observed median is ' num2str(diff_rtmedian) ' p = ' num2str(2*p_rtmedian)]);

% median accuracy
p_accmedian = length(find(bs_diffaccmedian>diff_accmedian))/nboots;
subplot(3,2,3);
bar(1:2,[magacc, notmagacc]);
hold on;
box off;
% confidence intervals
errorbar2(1:2,[magacc,notmagacc],[magaccci notmagaccci],1,'r');
set(gca,'XTickLabel',{'magnet syns','rest of pop'});
ylabel('accuracy for classification task');
% title();
subplot(3,2,4);
% calculate percentile less than observed value
p_accmedian = length(find(bs_diffaccmedian>diff_accmedian))/nboots
% add mark for observed median
[counts bincenters] = hist(bs_diffaccmedian,30);
% turn it into a probability distribution
bar(bincenters,counts/nboots,1);
hold on;
% yvalue of correct bin from histogram
yval=counts(find(min(bincenters-diff_accmedian)))/nboots;
% plot a line around that height of bin (will be small!)
plot(diff_accmedian,yval-(yval/2):yval+(yval/2),'r','Linewidth',8);
box off;
xlabel('bootstrapped difference in median accs for classification');
ylabel('probability of difference in medians'); 
title(['observed median is ' num2str(diff_accmedian) ' p = ' num2str(2*p_accmedian)]);

% color distance
p_cdistmedian = 1-length(find(bs_diffcdistmedian>diff_cdistmedian))/nboots
subplot(3,2,5);
bar(1:2,[magcdist, notmagcdist]);
hold on;
box off;
% confidence intervals
errorbar2(1:2,[magcdist,notmagcdist],[magcdistci notmagcdistci],1,'r');
set(gca,'XTickLabel',{'magnet syns','rest of pop'});
ylabel('accuracy for classification task');
% title();
subplot(3,2,6);
% calculate percentile less than observed value
p_cdistmedian = 1-length(find(bs_diffcdistmedian>diff_cdistmedian))/nboots
% add mark for observed median
[counts bincenters] = hist(bs_diffcdistmedian,30);
% turn it into a probability distribution
bar(bincenters,counts/nboots,1);
hold on;
% yvalue of correct bin from histogram
yval=counts(find(min(bincenters-diff_cdistmedian)))/nboots;
% plot a line around that height of bin (will be small!)
plot(diff_cdistmedian,yval-(yval/2):yval+(yval/2),'r','Linewidth',20);
box off;
xlabel('bootstrapped difference in median cdists for classification');
ylabel('probability of difference in medians'); 
title(['observed median is ' num2str(diff_cdistmedian) ' p = ' num2str(2*p_cdistmedian)]);



