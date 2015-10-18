%% s_ageHistograms
%
% Script to plot age histogram of synestehetes from Eagleman's database,
% grouping by those with many matches to fisher price toy versus those with
% fewer matches.

%% Navigate
%cd('~/projects/synesthesia/Eagleman')

%% load Eagleman RGB (small database)

% User ID and RGB colors come from this database of 6588 synesthetes with
% CG synesthesia
a = load('EaglemanColoredAlphabets');
uid = round(a.u_rgb(:,1));


%% Load Subject data (large database)

% Subject data (age, etc) come from database of all users (54450 users). We
% need to load up this data and find the correspondence between the two
% databases.

if ~exist('AboutPt', 'var')
    % only run it if the variables are not in the workspace because this
    % script  is slow to run
    run('LoadAboutPt')
end

% Note that AboutPtCols has one column for 'users_id' (col 1) and one for
% 'BatteryId' (col 2). Not sure which is right, that is which column of
% numbers should match up with uid from small database. This is obviosuly
% important because if the correspondence is wrong the analysis is
% meaningless.
%
% I am going with BatteryId because it contains every number in uid (from
% small database) whereas users_id is missing about half the numbers. The
% two columns contain similar but not identical numbers.

allusers.age = round(cell2mat(AboutPt(:,3)));
allusers.uid = round(cell2mat(AboutPt(:,2))); % assuming this right

allusers.uid2 = round(cell2mat(AboutPt(:,1)));  % assuming this is not right

% compare
% figure; plot(allusers.uid, allusers.uid2, '.')  

%% Match the two databases

% the function intersect sorts the data, so we have to do some somersaults
[dummy, ia, ib] = intersect(allusers.uid, uid);
[dummy, inds] = sort(ib);
ia = ia(inds);
age = allusers.age(ia);

% Check that we did it right
hasAge = find(isfinite(age));
assert(isequal(age(hasAge), allusers.age(ia(hasAge))))
assert(isequal(uid(hasAge), allusers.uid(ia(hasAge))))

%% Find good matches to toy

% rgb colors from database and toy and most common matches
rgb.eagle   = a.u_rgb(:,2:4, :);  % matrix of rgb values (n subj x rgb x 26 letters)
rgb.magnet  = fpSimulateData(size(rgb.eagle,1),'magnets');
rgb.common  = fpSimulateData(size(rgb.eagle,1),'most frequent');

rgb.magnet  = permute(rgb.magnet, [1 3 2]); % put the dimensionality in the same order as EAGLEMAN rgb (n, 3, 26)
rgb.common  = permute(rgb.common, [1 3 2]); % put the dimensionality in the same order as EAGLEMAN rgb (n, 3, 26)

% rgb to labels
labels.eagleman = fpRGB2ColorsJW(rgb.eagle);   % labels for Eagleman data
labels.magnet   = fpRGB2ColorsJW(rgb.magnet);  % labels for the magnet set
labels.common   = fpRGB2ColorsJW(rgb.common);  % labels for the common matches

% count  how many of EAGLEMAN RGBs agree with the magnet set and with most
% common matches
matches.magnet      = labels.eagleman == labels.magnet;
matches.common      = labels.eagleman == labels.common;
nummatches.magnet   = sum(matches.magnet, 2);
nummatches.common   = sum(matches.common, 2);

% find subjects with more than XX matches to toy / common matches
threshold =10;
thresh.magnet = nummatches.magnet > threshold;
thresh.common = nummatches.common > threshold;

%% Age histogram
bins = 15:5:50; % age bins
ages.magnet.matches    = hist(age(thresh.magnet), bins);
ages.magnet.nonmatches = hist(age(~thresh.magnet), bins);
ages.common.matches    = hist(age(thresh.common), bins);
ages.common.nonmatches = hist(age(~thresh.common), bins);

% --- Matches to magnet set --------------------------------------------
fH = figure(101); 
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 20);
[ax,h1,h2] = plotyy(2012-bins, ages.magnet.matches, ...
    2012-bins, ages.magnet.nonmatches);

set(ax(1), 'FontSize', 20, 'xlim', [1960 2000], 'XTick', 1930:10:2020);
set(ax(2), 'FontSize', 20, 'xlim', [1960 2000], 'XTick', 1930:10:2020);
set(h1, 'LineWidth', 2, 'Marker', 'o')
set(h2, 'LineWidth', 2, 'Marker', 'o')

legend({sprintf('%d or more matches', threshold+1),...
    sprintf('less than %d matches',  threshold+1)}, ...
    'Location', 'Best');
xlabel('Year of birth')
ylabel('Number')
title('Matches to toy')
saveas(fH, 'ageHistogramsMagnets.png')

% --- Matches to most common matches ----------------------------------
fH = figure(102); 
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 20);
[ax,h1,h2] = plotyy(2012-bins, ages.common.matches, ...
    2012-bins, ages.common.nonmatches);

set(ax(1), 'FontSize', 20, 'xlim', [1960 2000], 'XTick', 1930:10:2020);
set(ax(2), 'FontSize', 20, 'xlim', [1960 2000], 'XTick', 1930:10:2020);
set(h1, 'LineWidth', 2, 'Marker', 'o')
set(h2, 'LineWidth', 2, 'Marker', 'o')

legend({sprintf('%d or more matches', threshold+1),...
    sprintf('less than %d matches',  threshold+1)}, ...
    'Location', 'Best');
xlabel('Year of birth')
ylabel('Number')
title('Matches to most common matches')
saveas(fH, 'ageHistogramsCommon.png')


% --- Fraction of matches to magnet set and most common matches --------
% fH = figure(103); 
figure('name','matches to set and common by year','color',[1 1 1]);
set(gca, 'FontSize', 20);

fraction.magnet = ages.magnet.matches./...
    (ages.magnet.matches+ages.magnet.nonmatches);

fraction.common = ages.common.matches./...
    (ages.common.matches+ages.common.nonmatches);

plot(2010-bins, fraction.magnet);
% 
hold on;
% plot(2010-bins, fraction.common);

set(gca, 'FontSize', 20, 'xlim', [1960 2000], 'XTick', 1930:5:2020);


% title(sprintf('Fraction of subjects\nwith %d or more matches',
% threshold+1));
% legend({'Matches to magnet set' 'Matches to most common matches'}, ...
%     'Location', 'Best');
xlabel('Year of birth')
% legend boxoff;
box off;
plot2svg('matchesTOtemplatesbyyear.svg');


% saveas(fH, 'ageHistogramRelative.png')

%%
% for ii = 1:length(bins)-1
%     inds = age < bins(ii+1) & age >  bins(ii);
%     disp([bins(ii) sum(inds) median(nummatches(inds))])
%     medmatches(ii) = median(nummatches(inds));   
% end