% gc_population_model

% this is a model which takes a few bits of data and generates a prediction
% about what the population matching behavior of a group of grapheme color
% synesthetes will be like.

% let's assume that we ran EagleDBAnalysis73013.m

% let's grab indices to half the data set
pop_sample_indx = Shuffle(1:length(labels.eagleman));
% get test half
test_sample_indx = pop_sample_indx(end/2+1:end);
% keep first half as sample
pop_sample_indx =pop_sample_indx(1:end/2);

% let's see what our two halves look like
figure('Name','','Color',[1 1 1]);
subplot(1,2,1);
imagesc(permute(rgb.eagle(pop_sample_indx,:,:), [1 3 2]));

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title('sample data')

subplot(1,2,2);
imagesc(permute(rgb.eagle(test_sample_indx,:,:), [1 3 2]))
set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title('test data')


% now we need to find some things out about our population sample
% how many magnet synesthetes do we have?




% at the moment we would say that there are a set of contingencies between
% colors and letters.

% we estimate each one of these contingencies by just using the probability
% distribution of letters over colors.

% however different people are exposed to different contingencies.  broadly
% speaking we know of two
% 1. the contingencies everyone is exposed to in the culture
% 2. the magnet set

% so step 1 is determine the contingencies
% exposed to magnet set or not?

% if yes then for each letter flip a weighted coin where the weight is the
% probability that letter was assigned the magnet color

%     if it is not the magnet color, then flip a coin that is weighted with
%     by the probability distribution from the population

% if not in the magnet set then flip a coin weighted by the probability
% distribution from the population


% create a population of synesthetes and compare to original


% do we need competing models?  or could split half, so use half the data
% to build the model and the other half to test the prediction from the
% model.

