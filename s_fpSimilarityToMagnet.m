% Script. Gnerate many sets of RGB matches to the 26 letters, some of which
% are similar to the magnet set, some of which use the magnet set colors
% but assigned randomly to the letters, and some of which are purely random
% colors. Then calculate the similarity of these sets to the magnet set
% using Eagleman's distance metric.


% subset of subjects: magnet set colors with noise
m = fpSimulateData(1,'magnets');
RGB = fpSimulateData(50,'magnets noisy', 0.3);

% subset of subjects: shuffled colors (from magnet set)
RGB = cat(1, RGB, fpSimulateData(500,'shuffle'));

% subset of  subjects: random colors
RGB = cat(1, RGB, fpSimulateData(500,'random'));

%% visualize
figure(101);
image(RGB)
xlabel('letters')
ylabel('subjects')
title('Subjects')

figure(102);
image(m)
xlabel('letters')
title('Magnet set')

%% compute similarity
D = fpEaglemanDistance(m, RGB);
figure(103)
hist(D(:))