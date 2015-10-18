% Script to count the number of matches from Eagleman data set to general trends in the data
% . Compare this to the number of matches for shuffled data set, and
% plot.



% assumes you have or are running s_EagleDBAnalysis73013.m



%% what we want next is the most frequently occuring color for each letter.
%% for some letters this is essentially random, but for others it is pretty
%% systematic.

% make a matrix that has the subjects with exactly the most frequent
% matches
% we'll get as many rows as we have subjects, which makes comparison easier
rgb.fq = fpSimulateData(n,'most frequent');

% put the dimensionality in the same order as EAGLEMAN rgb (n, 3, 26)
% annoying as we unpermute it to get the image (haha).
rgb.fq = permute(rgb.fq, [1 3 2]);

%% RGB => Labels (0:10, based on NW's matrix)
labels.fq       = fpRGB2ColorsJW(rgb.fq);       % labels for the most frequent colors



%% Matches

% count up how many of EAGLEMAN RGBs agree with the frequently chosen
% colors
matches = labels.eagleman == labels.fq;
nummatches.eagle2fq  = sum(matches, 2);

% count up how many of Eagleman Shuffled by col match the most frequent
% colors
matches = labels.eagleShuffledByCol == labels.fq;
nummatches.eagleShuffleCol2fq = sum(matches,2);

% count up how many of Eagleman Shuffled by row match the most frequent
% colors
matches = labels.eagleShuffledByRow == labels.fq;
nummatches.eagleShuffleRow2fq = sum(matches,2);

% count up how many of Rich agree with most frequent
matches = labels.rich == labels.fq;
nummatches.rich2fq = sum(matches,2);

% already computed random

%% plot histograms

% matches to most frequent template
figure('name','number of matches to most frequent matches log scale', 'Color', [1 1 1],'Position',get(0,'ScreenSize'));

fq = hist(nummatches.eagle2fq, 0:26);
let = hist(nummatches.eagleShuffleCol2fq, 0:26);
sub = hist(nummatches.eagleShuffleRow2fq, 0:26);
ric = hist(nummatches.rich2fq, 0:26);
s(4,:) = hist(nummatches.random, 0:26);
uni = hist(nummatches.uniform, 0:26);


set(gcf, 'Color', 'w'); 
plot(0:26, fq, 'ro-', 0:26, sub, 'go-',0:26, let,'ko-', 0:26, ric, 'bo-',0:26,s(4,:),'co-',...
    0:26, uni, 'mo-','LineWidth', 2);
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Eagleman shuffled by subject','Eagleman shuffled by letter',...
    'rich dataset','random rgb','uniform labels');
legend boxoff;
title('Num matches to most frequent matches')
box off;
xlabel('number of matches to most frequent template');
ylabel('log number of synesthetes');

% matches to most frequent template% 
% % matches on a linear scale just for the heck of it
figure('name','number of matches to most frequent matches log scale', 'Color', [1 1 1]);

plot(0:26, fq, 'ro-', 0:26, sub, 'go-',0:26, let,'ko-', 0:26, ric, 'bo-',0:26,s(4,:),'co-',...
    0:26, uni, 'mo-','LineWidth', 2);set(gca, 'FontSize', 16, 'XLim', [0 26]);

legend('Eagleman RGB', 'Eagleman shuffled by subject','Eagleman shuffled by letter',...
    'rich dataset','random rgb','uniform labels');
legend boxoff;
title('Num matches to most frequent matches')
box off;
xlabel('number of matches to most frequent template');
ylabel('log number of synesthetes');


% 
% 
% % matches to real most frequent matches
% figure('name','matches to most fq minus shuffled matches to most fq', 'Color', [1 1 1]);
% 
% m = hist(nummatches.eagle2fq, 0:26);
% s = hist(nummatches.eagleShuffle2fq, 0:26);
% 
% set(gcf, 'Color', 'w'); 
% % plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
% set(gca, 'FontSize', 16, 'XLim', [0 26],'YScale','log');
% hold on;
% 
% 
% legend('Eagleman RGB')
% title('matches minus shuffled matches to most frequent matches')
% box off;
% plot(0:26,zeros(1,27),'k-','LineWidth',2);

% 
% 
% % matches to shuffled most frequent matches
% figure('name', 'matches to shuffled most frequent matches', 'Color', [1 1 1]);
% 
% m = hist(nummatches.eagle2fqShuffled, 0:26);
% s = hist(nummatches.eagleShuffle2fqShuffled, 0:26);
% 
% set(gcf, 'Color', 'w'); 
% plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
% set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');
% 
% legend('Eagleman RGB', 'Shuffled RGB')
% title('Num matches to shuffled most frequent matches')


%% plot colors for matches that are close to most frequent matches

% The best n matches (say, 500) ------------------------------------------
%   Rank the subjects by how many matches to letter set
% 
% % let's generate a bunch of these
% for i =1:7
% n = i*200;
% 
% [Y, ranking] = sort(nummatches.eagle2fq, 'descend');
% best = ranking(1:n);
% 
% figure('name', ['top ' num2str(n) ' matches to most frequent'], 'Color', [1 1 1]); 
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.fq(1:n/10,:,:); 
% theResult = rgb.eagle(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('Best matches to most frequent matches')
% end
% 
% % The best 500 matches to the shuffled most frequent matches -----------------------
% [Y, ranking] = sort(nummatches.eagle2fqShuffled, 'descend');
% best = ranking(1:n);
% 
% figure('name', ['top ' num2str(n) ' matches to shuffled most frequent'], 'Color', [1 1 1]); 
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.fqShuffled(1:n/10,:,:); 
% theResult = rgb.eagle(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('Best matches to shuffled most frequent matches')


% if you wanted to compare it to the earlier figure, you would shuffle the
% data instead of the most frequent matches and then sort.  the most frequent matches is just one
% throw instead of n?
% 
% % All the shuffled matches 
% figure('name', ['top ' num2str(n) ' shuffled matches to most frequent'], 'Color', [1 1 1]); 
% n = 500;
% [Y, ranking] = sort(nummatches.eagleShuffle2fq, 'descend');
% best = ranking(1:n);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.fq(1:round(n/10),:,:); 
% theResult = rgb.eagleShuffled(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('All shuffled matches to most frequent matches')
% 
% 
% 
% % All the matches
% figure(5); clf
% n = size(rgb.eagle,1);
% [Y, ranking] = sort(nummatches.eagle2fq, 'descend');
% best = ranking(1:n);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.fq(1:round(n/10),:,:); 
% theResult = rgb.eagle(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('All matches to most frequent matches')
% 
% % All the shuffled matches 
% figure(6); clf
% n = size(rgb.eagle,1);
% [Y, ranking] = sort(nummatches.eagleShuffle2fq, 'descend');
% best = ranking(1:n);
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.fq(1:round(n/10),:,:); 
% theResult = rgb.eagleShuffled(best,:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]))
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
% title('All shuffled matches to most frequent matches')
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
% 
