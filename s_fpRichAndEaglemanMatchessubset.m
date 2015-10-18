
% let's assume that you have already run s_fpRichAndEaglemanMatches.m

% let's just pull oot the subset of things we want from the output of that
% script.  this is a bit retarded as I could have done this in the first
% place....

rgbsubset.eagle = rgb.eagle(:,:,nooverlap);
rgbsubset.eagleShuffledByRow = rgb.eagleShuffledByRow(:,:,nooverlap);
rgbsubset.eagleShuffledByCol = rgb.eagleShuffledByCol(:,:,nooverlap);
rgbsubset.magnets = rgb.magnets(:,:,nooverlap);
rgbsubset.fq = rgb.fq(:,:,nooverlap);

labelssubset.eagleman = labels.eagleman(:,nooverlap);
labelssubset.eagleShuffledByRow = labels.eagleShuffledByRow(:,nooverlap);
labelssubset.eagleShuffledByCol = labels.eagleShuffledByCol(:,nooverlap);
labelssubset.magnets = labels.magnet(:,nooverlap);
labelssubset.rich = labels.rich(:,nooverlap);
labelssubset.fq = labels.fq(:,nooverlap);


richsubset.frequencies = rich.frequencies(nooverlap,:);
richsubset.cumFrequencies = rich.cumFrequencies(nooverlap,:);


%% Matches for magnets

% count up how many of EAGLEMAN RGBs agree with the magnet set
matches = labelssubset.eagleman == labelssubset.magnets;
nummatchessubset.eagle = sum(matches, 2);

% same for shuffled EAGLEMAN RGB (by rows)
matches  = labelssubset.eagleShuffledByRow == labelssubset.magnets;
nummatchessubset.eagleShuffleByRow  = sum(matches, 2);

% same for shuffled EAGLEMAN RGB (by cols)
matches  = labelssubset.eagleShuffledByCol == labelssubset.magnets;
nummatchessubset.eagleShuffleByCol  = sum(matches, 2);

% same for Rich distribution
matches  = labelssubset.rich == labelssubset.magnets;
nummatchessubset.rich = sum(matches, 2);

%% Matches for frequency

% count up how many of EAGLEMAN RGBs agree with the most frequent colors
matches = labelssubset.eagleman == labelssubset.fq;
nummatchessubset.eagle2fq = sum(matches, 2);

% same for shuffled EAGLEMAN RGB (by rows)
matches  = labelssubset.eagleShuffledByRow == labelssubset.fq;
nummatchessubset.eagleShuffleByRow2fq  = sum(matches, 2);

% same for shuffled EAGLEMAN RGB (by cols)
matches  = labelssubset.eagleShuffledByCol == labelssubset.fq;
nummatchessubset.eagleShuffleByCol2fq  = sum(matches, 2);

% same for Rich distribution
matches  = labelssubset.rich == labelssubset.fq;
nummatchessubset.rich2fq = sum(matches, 2);



%% plot histograms

% matches to most frequent template
figure('name','number of matches to most frequent matches log scale', 'Color', [1 1 1]);

fq = hist(nummatchessubset.eagle2fq, 0:length(no_o_letters));
let = hist(nummatchessubset.eagleShuffleByCol2fq, 0:length(no_o_letters));
sub = hist(nummatchessubset.eagleShuffleByRow2fq, 0:length(no_o_letters));
ric = hist(nummatchessubset.rich2fq, 0:length(no_o_letters));
% 
% disp('most frequent hist');
% [1:26' m]

set(gcf, 'Color', 'w'); 
plot(0:length(no_o_letters), fq, 'ro-', 0:length(no_o_letters), sub, 'go-',0:length(no_o_letters), let,'ko-', 0:length(no_o_letters), ric, 'bo-', 'LineWidth', 2);
% plot(0:length(no_o_letters), m-s, 'ro-', 'LineWidth', 2)
set(gca, 'FontSize', 16, 'XLim', [0 length(no_o_letters)], 'YScale', 'log');

legend('Eagleman RGB', 'Eagleman shuffled by subject','Eagleman shuffled by letter')
legend boxoff;
title('Num matches to most frequent matches')
box off;
xlabel('number of matches to most frequent template');
ylabel('log number of synesthetes');

% matches to most frequent template% 
% % matches on a linear scale just for the heck of it
figure('name','number of matches to most frequent matches log scale', 'Color', [1 1 1]);

plot(0:length(no_o_letters), fq, 'ro-', 0:length(no_o_letters), sub, 'go-',0:length(no_o_letters), let,'ko-', 0:length(no_o_letters), ric, 'bo-', 'LineWidth', 2);
set(gca, 'FontSize', 16, 'XLim', [0 length(no_o_letters)]);

legend('Eagleman RGB', 'Eagleman shuffled by subject','Eagleman shuffled by letter')
legend boxoff;
title('Num matches to most frequent matches')
box off;
xlabel('number of matches to most frequent template');
ylabel('log number of synesthetes');





% matches to real magnet set, eagleman v shuffled

m = hist(nummatchessubset.eagle, 0:length(no_o_letters));
sset(1,:) = hist(nummatchessubset.eagleShuffleByRow, 0:length(no_o_letters));
sset(2,:) = hist(nummatchessubset.eagleShuffleByCol, 0:length(no_o_letters));
sset(3,:) = hist(nummatchessubset.rich, 0:length(no_o_letters));

controlCond{1} = 'Eagleman RGB, shuffled within subject';
controlCond{2} = 'Eagleman RGB, shuffled within letter';
controlCond{3} = 'Rich distribution';
%
% for plotnum = 1:3
%     figure(plotnum);clf
%
%     set(gcf, 'Color', 'w');
%     plot(0:26, m, 'ro-', 0:26, s(plotnum,:), 'go-', 'LineWidth', 2)
%     %plot(0:26, m-s, 'ro-', 'LineWidth', 2)
%     set(gca, 'FontSize', 16, 'XLim', [0 26], 'YScale', 'log');
%
%     legend('Eagleman RGB', controlCond{plotnum})
%     title('Num matches to magnet set')
%
% end
%

% let's put these on a single plot
figure('name','compare magnet syns to rich distribution and shuffled eagleman',...
    'Color',[1 1 1]);

colors = [0 1 0; 0 .5 0; 0 0 1];

plot(0:length(no_o_letters), m, 'ro-', 0:length(no_o_letters), sset(1,:), 'go-',0:length(no_o_letters), sset(2,:), 'ko-',...
    0:length(no_o_letters), sset(3,:), 'bo-','LineWidth', 2);
%plot(0:length(no_o_letters), m-s, 'ro-', 'LineWidth', 2)

set(gca, 'FontSize', 16, 'XLim', [0 length(no_o_letters)], 'YScale', 'log');

title('Num matches to magnet set')
box off;
legend('Eagleman','Eagleman shuffled by subject',...
    'Eagleman shuffled by letter','Rich simulated dataset',...
    'Location','NorthEast')
legend boxoff;
ylabel('log number of subjects');
xlabel('number of matches to magnet set');



% let's put these on a single plot
figure('name','compare magnet syns to rich distribution and shuffled eagleman',...
    'Color',[1 1 1]);

colors = [0 1 0; 0 .5 0; 0 0 1];

plot(0:length(no_o_letters), m, 'ro-', 0:length(no_o_letters), sset(1,:), 'go-',0:length(no_o_letters), sset(2,:), 'ko-',...
    0:length(no_o_letters), sset(3,:), 'bo-','LineWidth', 2);

set(gca, 'FontSize', 16, 'XLim', [0 length(no_o_letters)]);

title('Num matches to magnet set')
box off;
legend('Eagleman','Eagleman shuffled by subject',...
    'Eagleman shuffled by letter','Rich simulated dataset',...
    'Location','NorthEast')
legend boxoff;
ylabel('log number of subjects');
xlabel('number of matches to magnet set');













% % then do the differential between the eagleman and rich datasets
% 
% % subtract two vectors
%  EvsR = (m-ss(3,:));
% %  find positive values (kind of hack, but these are to the right)
% msyns=sum(EvsR(find(EvsR>0)));
% 
% 
% 
% % let's put these on a single plot
% figure('name','matches to magnet set :eagleman minus rich simulated data',...
%     'Color',[1 1 1]);
% 
%     plot(0:length(nooverlap), m-ss(3,:), 'ro-', 'LineWidth', 2)
% 
% set(gca, 'FontSize', 16, 'XLim', [0 length(nooverlap)], 'YScale', 'log');
% 
% title('Num matches to magnet set')
% box off;
% ylabel('log number of subjects');
% xlabel('number of matches to magnet set');


