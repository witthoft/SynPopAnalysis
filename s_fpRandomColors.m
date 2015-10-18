% script to muck around with simluated Eagleman data

% number of subjects to simulate
n = 50;

RGB = []; txt = [];

% magnet set colors
RGB{1} = fpSimulateData(n,'magnets');
txt{1} = 'magnet set';

% magnet set colors with noise
RGB{2} = fpSimulateData(n,'magnets noisy', 0.1);
txt{2} = 'magnet set with noise (sd = 0.1)';

% magnet set colors with more noise
RGB{3} = fpSimulateData(n,'magnets noisy', 0.3);
txt{3} = 'magnet set with noise (sd = 0.3)';

% shuffled colors (from magnet set)
RGB{4} = fpSimulateData(n,'shuffle');
txt{4} = 'shuffled colors';

% shuffled with noise (from magnet set)
RGB{5} = fpSimulateData(n,'shuffle noisy', 0.1);
txt{5} = 'shuffled with noise (sd = 0.1)';

% random colors (from uniform distribution)
RGB{6} = fpSimulateData(n,'rand');
txt{6} = 'random colors';


n = numel(RGB);
nrows = round(sqrt(n));
ncols = ceil(n/nrows);

%% look at them

figure(101); 
for ii = 1:n
    subplot(nrows, ncols,ii)
    image(RGB{ii})
    xlabel('letters')
    ylabel('subjects')
    title(txt{ii})
end

%%
figure(102); 
for ii = 1:n
    subplot(nrows, ncols,ii)
    D = fpEaglemanDistance(RGB{1}, RGB{ii});
    hist(D(:), 0:.1:3)
    xlim([0 3])
    title(txt{ii})
end
