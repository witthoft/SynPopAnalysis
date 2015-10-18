%% s_BootstrapCommonMatches

% Which letters in Eagleman DB are matched to a color significantly more
% frequently than would be expected by chance?

%% load eagleman rgb
a = load('EaglemanColoredAlphabets');
rgb.eagle = a.u_rgb(:,2:4, :);  % matrix of rgb values (n subj x rgb x 26 letters)

%% clean up values outside the range [0 1]
rgb.eagle(rgb.eagle>1) = 1; rgb.eagle(rgb.eagle<0) = 0; 

letters = a.labels;     % a - z
colors = {'black' 'white' 'red' 'green' 'yellow' 'blue' 'brown' 'purple' 'pink' 'orange' 'grey'};
clear a;                % it's 11 MB, so save some memory
n = size(rgb.eagle, 1); % number of subjects


%% Get the most commonly matched colors. 

% we'll get as many rows as we have subjects, which makes comparison easier
rgb.common  = fpSimulateData(n,'common');

% put the dimensionality in the same order as EAGLEMAN rgb (n, 3, 26)
rgb.common  = permute(rgb.common, [1 3 2]);

%% RGB => Labels (0:10, based on NW's matrix)
labels.eagleman = fpRGB2ColorsJW(rgb.eagle);              % labels for Eagleman data
labels.common   = fpRGB2ColorsJW(rgb.common);            % labels for the magnet set

%% Count how often each color is matched to each letter

% Make a 26 by 11 matrix. Each cell will indicate the number of times 
% letter i was matched to color j
colorFrequenciesByLetter = zeros(26,11);

% loop acorss the 26 letters
for i = 1:26 
    % loop across the 11 colors
    for j = 0:10

        % count the number of times letter 'i' was matched to color 'j'
        colorFrequenciesByLetter(i,j+1) = ...
            sum(labels.eagleman(:,i) == j);
    end
end

% divide by num subjects to convert to probability ([0 1])
colorFrequenciesByLetter = colorFrequenciesByLetter/n;

% get the highest frequency color match for each letter (e.g, freq A=>red)
[maxColorFrequenciesByLetter, whichColor] = max(colorFrequenciesByLetter, [], 2);

%% Simulate matches based on frequencies of colors in the whole data set

% count how often each color appears in the whole dataset
colorFrequencies = zeros(1,11);
for c = 0:10
    colorFrequencies(c+1) = sum(labels.eagleman(:) == c)/n/26;
end
cummulativeFrequencies = cumsum(colorFrequencies);

% now make nboot alphabets based on these 26 frequencies 
nboot      = 1000;
tmp        = rand(n ,nboot); % random numbers which will be converted to labels
labels.sim = zeros(n, nboot, 'uint8'); % this will hold the labels

% trick for getting the labels right
for c = 1:10
    labels.sim(tmp > cummulativeFrequencies(c)) = c;
end

% count how often each color appears in the boostrapped distribution
colorFrequenciesSim = zeros(11, nboot);
for ii = 1:nboot
    for c = 0:10
        colorFrequenciesSim(c+1,ii) = sum(labels.sim(:,ii) == c)/n;
    end
end

maxColorSim = max(colorFrequenciesSim,[], 1);
prcnts = prctile(maxColorSim, [2.5, 50, 97.5])

%%
figure;
bins = 0:.05:1;
m1 = hist(maxColorSim, bins);
m2 = hist(maxColorFrequenciesByLetter, bins);
figure(1); plot(bins, m1/sum(m1)*26, bins, m2/sum(m2)*26, 'o-', 'LineWidth', 2)


