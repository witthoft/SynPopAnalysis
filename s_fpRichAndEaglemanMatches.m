% Script to count the number of matches from Eagleman data set and magnet
% set. Compare this to the number of matches for shuffled data set, and
% plot.

saveFigures = false;


%% load up rich data
tmp = load('RichMatches');
% this loads the variable rich, which is a matrix 26x11 (letters x
% frequency of color matches)
% though the rows sum very close to 100, they don't quite usually.  mostly
% the error is less than one, but there are some exceptions, not sure why.
% maybe some matches couldn't be assigned to any of the color categories?

% the colors are stored in a different order than our labels
% from their paper
rich2eagleman = [1 3 6 8 4 9 11 10 7 5 2];

% letter x color matrix, in same order as our color labels
rich.frequencies = tmp.rich(:, rich2eagleman);



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
rgb.eagleShuffledByRow = rgb.eagle;
rgb.eagleShuffledByCol = rgb.eagle;
% % make a set of random rgb values
rgb.random = rand(size(rgb.eagle));

% shuffle by column
% loop over number of subjects and shuffle each one's letter-color mappings

for ii = 1:26
    rgb.eagleShuffledByCol(:,:,ii) = rgb.eagle(Shuffle(1:n), :, ii);
end

% We can shuffle by row (across subjects, within letters) or column (across
% letters, within subject)

% shuffle by row
% loop over number of letters and shuffle each matches

for ii = 1:n
    % Shuffle all letters except A. A is usually red, so the null
    % distribution should maybe include that.
    % %     shuffleThese = setdiff(1:26, 1);
    shuffleThese = 1:26;
    rgb.eagleShuffledByRow(ii,:,shuffleThese) = rgb.eagle(ii, :, Shuffle(shuffleThese));
end

%% Get the pre-defined magnet set colors.

% we'll get as many rows as we have subjects, which makes comparison easier
rgb.magnets = fpSimulateData(n,'magnets');

% put the dimensionality in the same order as EAGLEMAN rgb (n, 3, 26)
rgb.magnets = permute(rgb.magnets, [1 3 2]);

%% RGB => Labels (0:10, based on NW's matrix)
labels.eagleman            = fpRGB2ColorsJW(rgb.eagle);              % labels for Eagleman data
labels.eagleShuffledByRow  = fpRGB2ColorsJW(rgb.eagleShuffledByRow); % labels for Eagleman data shuffled
labels.eagleShuffledByCol  = fpRGB2ColorsJW(rgb.eagleShuffledByCol); % labels for Eagleman data shuffled
labels.magnet              = fpRGB2ColorsJW(rgb.magnets);            % labels for the magnet set
labels.random              = fpRGB2ColorsJW(rgb.random);             % labels for randomly generated rgbs
% % how about a truly random data set with all rgbs randomly generated
% % however here the distribution of labels isn't uniform but proportional
% to how much of the rgb space each label occupies
% 
% figure('Color',[1 1 1]);
% hist(labels.random(:),[0:10])
% box off
% xlabel('color labels');
% ylabel('count');
% % %  we could also make one that is uniform in the label space
labels.uniform = randi(11,size(labels.eagleman))-1;
% figure('Color',[1 1 1]);
% hist(labels.uniform(:),[0:10])
% box off
% xlabel('color labels');
% ylabel('count');
% % have labels.uniform alreadfy
% 
% 
% 

%% Create a null distribution using rich data
tmp                   = rand(n , 26); % random numbers which will be converted to labels
labels.rich           = zeros(n, 26); % this will hold the labels
rich.cumFrequencies   = cumsum(rich.frequencies, 2);
rich.cumFrequencies   = bsxfun(@rdivide, rich.cumFrequencies, sum(rich.frequencies,2));

for ll = 1:26
    for ii = 1:n
        labels.rich(ii,ll) = find(rich.cumFrequencies(ll,:) >= tmp(ii,ll),1) - 1;
    end
end

%% Matches

% count up how many of EAGLEMAN RGBs agree with the magnet set
matches = labels.eagleman == labels.magnet;
nummatches.eagle = sum(matches, 2);

% same for shuffled EAGLEMAN RGB (by rows)
matches  = labels.eagleShuffledByRow == labels.magnet;
nummatches.eagleShuffleByRow  = sum(matches, 2);

% same for shuffled EAGLEMAN RGB (by cols)
matches  = labels.eagleShuffledByCol == labels.magnet;
nummatches.eagleShuffleByCol  = sum(matches, 2);

% same for Rich distribution
matches  = labels.rich == labels.magnet;
nummatches.rich = sum(matches, 2);

% same for randomly generated rgb data
matches = labels.random == labels.magnet;
nummatches.random = sum(matches, 2);

% same for uniform distribution of color labels (equivalent to binomial
% with k = 11)
matches = labels.uniform == labels.magnet;
nummatches.uniform = sum(matches, 2);



%% plot histograms

% matches to real magnet set, eagleman v shuffled

m = hist(nummatches.eagle, 0:26);
s(1,:) = hist(nummatches.eagleShuffleByRow, 0:26);
s(2,:) = hist(nummatches.eagleShuffleByCol, 0:26);
s(3,:) = hist(nummatches.rich, 0:26);
s(4,:) = hist(nummatches.random, 0:26);
s(5,:) = hist(nummatches.uniform, 0:26);

Cond{1} = 'Eagleman RGB';
Cond{2} = 'Eagleman RGB, shuffled within subject';
Cond{3} = 'Eagleman RGB, shuffled within letter';
Cond{4} = 'Rich distribution';
Cond{5} = 'Random distribution of RGB';
Cond{6} = 'uniform distribution of labels';

% let's put these on a single plot with log y axis
figure('name','compare magnet syns to rich distribution and shuffled eagleman',...
    'Color',[1 1 1],'Position',get(0,'ScreenSize'));

colors = [0 1 0; 0 .5 0; 0 0 1];

plot(0:26, m, 'ro-', 0:26, s(1,:), 'go-',0:26, s(2,:), 'ko-',...
    0:26, s(3,:), 'bo-',0:26,s(4,:),'co-',0:26, s(5,:), 'mo-',...
    'LineWidth', 2)
%plot(0:26, m-s, 'ro-', 'LineWidth', 2)

set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');

title('Num matches to magnet set')
box off;
legend(Cond,'Location','NorthEast')
legend boxoff;
ylabel('log number of subjects');
xlabel('number of matches to magnet set');

% 


% make a nice svg plot
%  plot2svg('magnetsVSnulls.svg');






% 
% % let's put these on a single plot with linear y axis
% figure('name','compare magnet syns to rich distribution and shuffled eagleman',...
%     'Color',[1 1 1]);
% 
% colors = [0 1 0; 0 .5 0; 0 0 1];
% 
% 
% plot(0:26, m, 'ro-', 0:26, s(1,:), 'go-',0:26, s(2,:), 'ko-',...
%     0:26, s(3,:), 'bo-',0:26,s(4,:),'co-',0:26, s(5,:), 'mo-',...
%     'LineWidth', 2)
% 
% set(gca, 'FontSize', 16, 'XLim', [0 26]);
% 
% title('Num matches to magnet set')
% box off;
% legend(Cond,'Location','NorthEast')
% 
% legend boxoff;
% ylabel('log number of subjects');
% xlabel('number of matches to magnet set');

% 
% % then do the differential between the eagleman and rich datasets
% 
% % subtract two vectors
% EvsR = (m-s(3,:));
% %  find positive values (kind of hack, but these are to the right)
% msyns=sum(EvsR(find(EvsR>0)));
% 
% 
% 
% % let's put these on a single plot
% figure('name','matches to magnet set :eagleman minus rich simulated data',...
%     'Color',[1 1 1]);
% 
% plot(0:26, m-s(3,:), 'ro-', 'LineWidth', 2)
% 
% set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');
% 
% title([num2str(msyns) '  more magnet syns than expected from Rich Data']);
% box off;
% ylabel('log number of subjects');
% xlabel('number of matches to magnet set');
% 
% % then do the differential between the eagleman and column shuffled
% 
% % subtract two vectors
% EvsCS = (m-s(2,:));
% %  find positive values (kind of hack, but these are to the right)
% msyns=sum(EvsCS(find(EvsCS>0)));
% 
% 
% 
% % let's put these on a single plot
% figure('name','  matches to magnet set :eagleman minus column shuffled data',...
%     'Color',[1 1 1]);
% 
% plot(0:26, m-s(2,:), 'ro-', 'LineWidth', 2)
% 
% set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');
% 
% title([num2str(msyns) 'more magnet syns than expected from shuffled columns']);
% box off;
% ylabel('log number of subjects');
% xlabel('number of matches to magnet set');

% 
% 
% 
% 
% if saveFigures,
%     thedir = '/Volumes/winawer/projects/synesthesia/Eagleman';
%     saveas(1, fullfile(thedir, 'histoEaglemanMatches'), 'fig')
%     saveas(1, fullfile(thedir, 'histoEaglemanMatches'), 'png')
%     
%     %     saveas(2, fullfile(thedir, 'histoShuffledMatches'), 'fig')
%     %     saveas(2, fullfile(thedir, 'histoShuffledMatches'), 'png')
%     %
%     figure(1); set(gca, 'YScale', 'linear', 'YLim', [0 100])
%     saveas(1, fullfile(thedir, 'histoLinEaglemanMatches'), 'fig')
%     saveas(1, fullfile(thedir, 'histoLinEaglemanMatches'), 'png')
%     
%     %     figure(2); set(gca, 'YScale', 'linear', 'YLim', [0 100])
%     %     saveas(2, fullfile(thedir, 'histoLinShuffledMatches'), 'fig')
%     %     saveas(2, fullfile(thedir, 'histoLinShuffledMatches'), 'png')
% end
% 
% 


% want to do a correlation between the frequencies in the rich data and the
% eagleman data
% rich frequencies with same order as eagleman
richfq = (rich.frequencies/100)';
% eagleman frequencies
% for each letter
eaglefq = [];
for i=1:26
%     get the number of matches to each color category
    [counts, bins]=hist(dbNumbered(:,i),11);
%     turn into a percent
   eaglefq(:,i)  = counts/length(dbLabeled);
    
end


% now we can unpack the matrices and do a straight correlation.  could plot
% letters with every letter colored.

richfqvect = richfq(:)';
eaglefqvect = eaglefq(:)';
% should do the same with letter matrix if we want it

figure('Name','Eagleman vs Rich Letter-Color Matching frequencies','Color',[1 1 1],'Position',get(0,'ScreenSize'));
scatter(richfqvect,eaglefqvect,'r.');



[rho, p] = corr(richfqvect', eaglefqvect')













