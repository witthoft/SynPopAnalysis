% Script to count the number of matches from Eagleman data set and magnet
% set. Compare this to the number of matches for shuffled data set, and
% plot.

saveFigures = false;


cutoff=300

%% load up eagleman rgb
a = load('EaglemanColoredAlphabets');
rgb.eagle = a.u_rgb(:,2:4, :);  % matrix of rgb values (n subj x rgb x 26 letters)


%%
% clean up values outside the range [0 1]
rgb.eagle(rgb.eagle>1) = 1; rgb.eagle(rgb.eagle<0) = 0;

letters = a.labels;     % a - z
clear a;                % it's 11 MB, so save some memory
n = size(rgb.eagle, 1); % number of subjects

%% make a shuffled set for comparison
rgb.eagleShuffled = rgb.eagle;

% each subject's matches are shuffled, which is the same as shuffling the
% rows in a matrix which is subject's by matches
% loop over number of subjects and shuffle each one's letter-color mappings
for ii = 1:n
    
    % Shuffle all letters except A. A is usually red, so the null
    % distribution should maybe include that.
    shuffleThese = setdiff(1:26, 1);
    rgb.eagleShuffled(ii,:,shuffleThese) = rgb.eagle(ii, :, Shuffle(shuffleThese));
end

%% Get the pre-defined magnet set colors.

% we'll get as many rows as we have subjects, which makes comparison easier
% makes a matrix where every row is a match to the magnet set
rgb.magnets = fpSimulateData(n,'magnets');

% put the dimensionality in the same order as EAGLEMAN rgb (n, 3, 26)
rgb.magnets = permute(rgb.magnets, [1 3 2]);

% shuffle magnet set colors to see how well real data matches the shuffled
% toy
rgb.magnetsShuffled = rgb.magnets(:, :, Shuffle(1:26));

%% RGB => Labels (0:10, based on NW's matrix)
labels.eagleman       = fpRGB2ColorsJW(rgb.eagle);         % labels for Eagleman data
labels.eagleShuffled  = fpRGB2ColorsJW(rgb.eagleShuffled); % labels for Eagleman data shuffled
labels.magnet         = fpRGB2ColorsJW(rgb.magnets);       % labels for the magnet set
labels.magnetShuffled = fpRGB2ColorsJW(rgb.magnetsShuffled);

%% Matches

% count up how many of EAGLEMAN RGBs agree with the magnet set
matches = labels.eagleman == labels.magnet;
nummatches.eagle2magnet  = sum(matches, 2);

% same for shuffled EAGLEMAN RGB
matches  = labels.eagleShuffled == labels.magnet;
nummatches.eagleShuffle2magnet  = sum(matches, 2);

% count up how many of EAGLEMAN RGBs agree with the SHUFFLED magnet set
matches = labels.eagleman == labels.magnetShuffled;
nummatches.eagle2magnetShuffled = sum(matches, 2);

% same for shuffled EAGLEMAN RGB
matches = labels.eagleShuffled == labels.magnetShuffled;
nummatches.eagleShuffle2magnetShuffled = sum(matches, 2);

%% plot histograms

% matches to real magnet set
figure('name','number of matches to magnet set log scale', 'Color', [1 1 1]);

m = hist(nummatches.eagle2magnet, 0:26);
s = hist(nummatches.eagleShuffle2magnet, 0:26);

set(gcf, 'Color', 'w');
plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB')
title('Num matches to magnet set')
box off;


% matches to real magnet set
figure('name','matches to set minus shuffled matches to set', 'Color', [1 1 1]);

m = hist(nummatches.eagle2magnet, 0:26);
s = hist(nummatches.eagleShuffle2magnet, 0:26);

set(gcf, 'Color', 'w');
% plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26],'Ylim', [-340, 500],'YScale','log');
hold on;


legend('Eagleman RGB')
box off;
plot(0:26,zeros(1,27),'k-','LineWidth',2);


% let's use this as one way to set the limit on the magnet synesthetes
% so this is the difference between the number people with m-matches to the
% magnet set and shuffled data with m-matches to the magnet set
mVSs = m-s;

numMagSynsdiff = sum(mVSs(find(mVSs>0)))
title([ num2str(numMagSynsdiff) ' magnet syns']);




% matches to shuffled magnet set
figure('name', 'matches to shuffled magnet set', 'Color', [1 1 1]);

m = hist(nummatches.eagle2magnetShuffled, 0:26);
s = hist(nummatches.eagleShuffle2magnetShuffled, 0:26);

set(gcf, 'Color', 'w');
plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB')
title('Num matches to shuffled magnet set')


if saveFigures,
    thedir = '/Volumes/winawer/projects/synesthesia/Eagleman';
    saveas(1, fullfile(thedir, 'histoEaglemanMatches'), 'fig')
    saveas(1, fullfile(thedir, 'histoEaglemanMatches'), 'png')
    
    saveas(2, fullfile(thedir, 'histoShuffledMatches'), 'fig')
    saveas(2, fullfile(thedir, 'histoShuffledMatches'), 'png')
    
    figure(1); set(gca, 'YScale', 'linear', 'YLim', [0 100])
    saveas(1, fullfile(thedir, 'histoLinEaglemanMatches'), 'fig')
    saveas(1, fullfile(thedir, 'histoLinEaglemanMatches'), 'png')
    
    figure(2); set(gca, 'YScale', 'linear', 'YLim', [0 100])
    saveas(2, fullfile(thedir, 'histoLinShuffledMatches'), 'fig')
    saveas(2, fullfile(thedir, 'histoLinShuffledMatches'), 'png')
end
box off;

%% plot colors for matches that are close to magnet set
% so this a different way in main script

% 
% 
% [Y, ranking] = sort(nummatches.eagle2magnet, 'descend');
% best = ranking(1:cutoff);
% 
% figure('name', ['top ' num2str(cutoff) ' matches to letter set'], 'Color', [1 1 1]);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.magnets(1:cutoff/10,:,:);
% theResult = rgb.eagle(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('Best matches to magnet set')
% 


% there are a lot of ways of choosing n for the final figure



[Y, ranking] = sort(nummatches.eagle2magnetShuffled, 'descend');
best = ranking(1:cutoff);

% 
% % these figures have been moved
% figure('name', ['top ' num2str(cutoff) ' matches to shuffled letter set'], 'Color', [1 1 1]);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.magnetsShuffled(1:cutoff/10,:,:);
% theResult = rgb.eagle(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('Best matches to shuffled magnet set')


% if you wanted to compare it to the earlier figure, you would shuffle the
% data instead of the magnet set and then sort.  the magnet set is just one
% throw instead of n?

% % All the shuffled matches
% figure('name', ['top ' num2str(n) ' shuffled matches to letter set'], 'Color', [1 1 1]);
% 
% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% best = ranking(1:n);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.magnets(1:round(n/10),:,:);
% theResult = rgb.eagleShuffled(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('All shuffled matches to magnet set')

% 
% 
% % All the matches
% figure('name','all matches to magnet set','Color',[1 1 1]);
% n = size(rgb.eagle,1);
% [Y, ranking] = sort(nummatches.eagle2magnet, 'descend');
% best = ranking(1:n);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.magnets(1:round(n/10),:,:);
% theResult = rgb.eagle(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('All matches to magnet set')
% 
% % All the shuffled matches
% figure('name','all suffled matches to magnet set','Color',[1 1 1]);
% n = size(rgb.eagle,1);
% [Y, ranking] = sort(nummatches.eagleShuffle2magnet, 'descend');
% best = ranking(1:n);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.magnets(1:round(n/10),:,:);
% theResult = rgb.eagleShuffled(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('All shuffled matches to magnet set')
% 
% if saveFigures,
%     thedir = '/Volumes/winawer/projects/synesthesia/Eagleman';
%     saveas(3, fullfile(thedir, 'rgbMatchesToMagnet'), 'fig')
%     saveas(3, fullfile(thedir, 'rgbMatchesToMagnet'), 'png')
%     
%     saveas(4, fullfile(thedir, 'rgbMatchesToShuffle'), 'fig')
%     saveas(4, fullfile(thedir, 'rgbMatchesToShuffle'), 'png')
%     
%     saveas(5, fullfile(thedir, 'rgbAllMatches'), 'fig')
%     saveas(5, fullfile(thedir, 'rgbAllMatches'), 'png')
%     
%     saveas(6, fullfile(thedir, 'rgbAllShuffledMatches'), 'fig')
%     saveas(6, fullfile(thedir, 'rgbAllShuffledMatches'), 'png')
% end



% one other method of shuffling I wanted to try.
% basically create a large number of random templates.  this would be best
% done without replacement, but the probability of generating duplicates is
% pretty small.

numrandtemplates = 10000;

% put the output here
% is 27 because there can be 0 to 26 matches
randtempmatchhists = zeros(numrandtemplates,27);

for i=1:numrandtemplates
    %     could be compressed but for ease of reading
    % make a template
    randtemp= round(10*rand(1,26));
    %     match to size of data set
    randtemp=repmat(randtemp(1,:),6588,1);
    %     check for number of matches between each subject and this template
    mtchs = labels.eagleman==randtemp;
%     get the histogram from the full matrix (i.e. the number of subjects
%     with n matches where n ranges from 0 to 26
    randtempmatchhists(i,:)=hist(sum(mtchs,2),0:26);
end
    
    
    
%  could plot the whole thing 
% figure;plot(0:26,randtempmatchhists)

% but lets do the mean and standard devation of each column

meanofrtm = mean(randtempmatchhists);
stdofrtm = std(randtempmatchhists);

% this is kind of sneaking up on the binomial distribution approach from
% the data direction.


    
% then various other ways of comparing the data when shuffled and not.
% these figures are nice, but we still lack a way of drawing a boundary
% between those who had the set and those who didn't.  we could say 10 and
% that would be pretty safe.  or we could say that given that forgetting a
% and c, how many matches out of 24 should we see?

% maybe we should superimpose the other kinds of null distribution on the
% data set?

% plot histograms


% make simulated data set

simdata = randi(11,size(matches))-1;
% should make a simulated data figure...
% count up matches to magnet set
% the labels variable comes from s_fpEaglemanMatches
nsimmatches = simdata == labels.magnet;
% then sum across the columns
simletmatches = sum(nsimmatches,2);
% then add to the plot jon made


% matches to real magnet set
figure('name','number of matches to magnet set log scale', 'Color', [1 1 1]);

m = hist(nummatches.eagle2magnet, 0:26);
s = hist(nummatches.eagleShuffle2magnet, 0:26);
k = hist(simletmatches, 0:26);

set(gcf, 'Color', 'w');
plot(0:26, m, 'ro-', 0:26, s, 'go-', 0:26, k, 'bo-','LineWidth', 2)
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB','simulated data')
title('Num matches to magnet set')
box off;

% as expected, the simulated data curve is shifted to the left of the
% shuffled data curve, because the data themselves are not random in the
% way that we imagined for the psychscience paper.

% suppose we fix the first column to be red say 50 % of the time, and then
% make c yellow about 40 % of the time.


% so for the first column we make .4*6588 ~=2635  length vector of 2s.
a= 2*ones(1,2635);
% the rest of the vector can be uniformly sampled from the other
% categories.  in this case the range can be 1-10 and then all the 2s set
% to 0.
b=randi(10,[1,6588-2635]);
% set 2s = 0
b(find(b==2))=0;
% combine the vectors
acol = shuffle([a b]);
%  figure;hist(acol,11) looks good
% lets do the same for c
c = 4*ones(1,1976);
crest = randi(10,[1,6588-1976]);
c(find(c==4))=0;
ccol = shuffle([c crest]);
% now lets add them to our simulated data

adjsimdata = simdata;
adjsimdata(:,1)=acol';
adjsimdata(:,3)=ccol';

% now we can see how well this fits the shuffled data we actually have
% count up matches to magnet set
% the labels variable comes from s_fpEaglemanMatches
adjsimmatches = adjsimdata == labels.magnet;
% then sum across the columns
adjsimletmatches = sum(adjsimmatches,2);
% then add to the plot jon made


% matches to real magnet set
figure('name','number of matches to magnet set log scale', 'Color', [1 1 1]);

m = hist(nummatches.eagle2magnet, 0:26);
s = hist(nummatches.eagleShuffle2magnet, 0:26);
k = hist(adjsimletmatches, 0:26);

set(gcf, 'Color', 'w');
plot(0:26, m, 'ro-', 0:26, s, 'go-', 0:26, k, 'bo-','LineWidth', 2)
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB','adjusted (a,c) simulated data')
title('Num matches to magnet set')
box off;

% this shifts it a little closer.  anyway, in the extreme we are just
% simulating the original dataset.  what if we just don't consider the
% first column (the a's) at all?









% matches to real magnet set
figure('name','number of matches to magnet set log scale', 'Color', [1 1 1]);

m = hist(nummatches.eagle2magnet, 0:26);
s = hist(nummatches.eagleShuffle2magnet, 0:26);
k = hist(simletmatches, 0:26);

set(gcf, 'Color', 'w');
plot(0:26, m, 'ro-', 0:26, s, 'go-', 0:26, k, 'bo-','LineWidth', 2)
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB','simulated data')
title('Num matches to magnet set')
box off;


