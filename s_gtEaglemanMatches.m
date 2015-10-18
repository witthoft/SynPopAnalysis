% Script to count the number of matches from Eagleman data set to general trends in the data
% . Compare this to the number of matches for shuffled data set, and
% plot.

saveFigures = false;



%% load up eagleman rgb
a = load('EaglemanColoredAlphabets');
rgb.eagle = a.u_rgb(:,2:4, :);  % matrix of rgb values (n subj x rgb x 26 letters)


%%
% clean up values outside the range [0 1]
rgb.eagle(rgb.eagle>1) = 1; rgb.eagle(rgb.eagle<0) = 0; 

letters = a.labels;     % a - z
clear a;                % it's 11 MB, so save some memory
n = size(rgb.eagle, 1); % number of subjects

% already have this from fp_visualizeShuffledEagDB.m

%% make a shuffled set for comparison 
% rgb.eagleShuffled = rgb.eagle; 
% 
% % loop over number of subjects and shuffle each one's letter-color mappings
% for ii = 1:n
%     
%     % Shuffle all letters except A. A is usually red, so the null
%     % distribution should maybe include that.
%     shuffleThese = setdiff(1:26, 1); 
%     rgb.eagleShuffled(ii,:,shuffleThese) = rgb.eagle(ii, :, Shuffle(shuffleThese));
% end


% so, some processes are duplicated in the various scripts.  there are
% already shuffled data in rgb from fp_visualizeShuffledEagDB.m so it is
% unneccesary to do it again
% rgb = 
% 
%                  eagle: [6588x3x26 double]
%     eagleShuffledByRow: [6588x3x26 double]
%     eagleShuffledByCol: [6588x3x26 double]
%                magnets: [6588x3x26 double]


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



% shuffle most frequent colors to see how well real data matches the shuffled
% frequent data
% now every subject is the same but the color matches are shuffled.
rgb.fqShuffled = rgb.fq(:, :, Shuffle(1:26));

% % this also might already exist
% %% RGB => Labels (0:10, based on NW's matrix)
% labels.eagleman       = fpRGB2ColorsJW(rgb.eagle);         % labels for Eagleman data
% labels.eagleShuffled  = fpRGB2ColorsJW(rgb.eagleShuffled); % labels for Eagleman data shuffled
labels.fq       = fpRGB2ColorsJW(rgb.fq);       % labels for the most frequent colors
labels.fqShuffled = fpRGB2ColorsJW(rgb.fqShuffled); 
% some other labels also already exist

% labels = 
% 
%               eagleman: [6588x26 double]
%     eagleShuffledByRow: [6588x26 double]
%     eagleShuffledByCol: [6588x26 double]
%                 magnet: [6588x26 double]
%                   rich: [6588x26 double]



%% Matches

% count up how many of EAGLEMAN RGBs agree with the frequently chosen
% colors
matches = labels.eagleman == labels.fq;
nummatches.eagle2fq  = sum(matches, 2);

% % same for shuffled EAGLEMAN RGB
matches  = labels.fqShuffled == labels.fq;
nummatches.fqShuffle2fq  = sum(matches, 2);
% 
% % count up how many of EAGLEMAN RGBs agree with the SHUFFLED frequent colors set
% matches = labels.eagleman == labels.fqShuffled;
% nummatches.eagle2fqShuffled = sum(matches, 2);
% 
% % same for shuffled EAGLEMAN RGB
% matches = labels.eagleShuffled == labels.fqShuffled;
% nummatches.eagleShuffle2fqShuffled = sum(matches, 2);

% what we really want is the column shuffled data
% matches = labels.


%% plot histograms

% matches to real magnet set
figure('name','number of matches to most frequent matches log scale', 'Color', [1 1 1]);

m = hist(nummatches.eagle2fq, 0:26);
s = hist(nummatches.fqShuffle2fq, 0:26);

disp('most frequent hist');
[1:26' m]

set(gcf, 'Color', 'w'); 
plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
% plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB')
title('Num matches to most frequent matches')
box off;

% 
% % matches on a linear scale just for the heck of it
% figure('name','number of matches to most frequent matches linear scale', 'Color', [1 1 1]);
% 
% m = hist(nummatches.eagle2fq, 0:26);
% s = hist(nummatches.eagleShuffle2fq, 0:26);
% 
% set(gcf, 'Color', 'w'); 
% plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
% % plot(0:26, m-s, 'ro-', 'LineWidth', 2)
% set(gca, 'FontSize', 16, 'XLim', [0 26], 'Ylim', [0 80]);
% 
% legend('Eagleman RGB', 'Shuffled RGB')
% title('Num matches to most frequent matches')
% box off;

% matches to real most frequent matches
figure('name','matches to most fq minus shuffled matches to most fq', 'Color', [1 1 1]);

% m = hist(nummatches.eagle2fq, 0:26);
% s = hist(nummatches.eagleShuffle2fq, 0:26);

set(gcf, 'Color', 'w'); 
% plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
plot(0:26, m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26],'YScale','log');
hold on;


legend('Eagleman RGB')
title('matches minus shuffled matches to most frequent matches')
box off;
plot(0:26,zeros(1,27),'k-','LineWidth',2);



% matches to shuffled most frequent matches
figure('name', 'matches to shuffled most frequent matches', 'Color', [1 1 1]);

% m = hist(nummatches.eagle2fqShuffled, 0:26);
% s = hist(nummatches.eagleShuffle2fqShuffled, 0:26);

set(gcf, 'Color', 'w'); 
plot(0:26, m, 'ro-', 0:26, s, 'go-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

legend('Eagleman RGB', 'Shuffled RGB')
title('Num matches to shuffled most frequent matches')


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
