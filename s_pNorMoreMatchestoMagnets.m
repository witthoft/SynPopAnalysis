% Script to count the number of matches from Eagleman data set and magnet
% set. Compare this to the number of matches for shuffled data set, and
% plot.
% have

% labels =
%
%               eagleman: [6588x26 double]
%     eagleShuffledByRow: [6588x26 double]
%     eagleShuffledByCol: [6588x26 double]
%                 magnet: [6588x26 double]
%                 random: [6588x26 double]
%                uniform: [6588x26 double]
%                   rich: [6588x26 double]




% do monte carlo to find out probability of observing n or more matches
% based on several choices of null distribution

sim.colshuffled = zeros(1,27);
sim.rowshuffled = zeros(1,27);
sim.rich = zeros(1,27);
sim.uniform = zeros(1,27);
sim.randomRGB = zeros(1,27);


numsims = 10;

% so for s shuffles
for i=1:numsims
    
    
    % column shuffled data
    % make a simulated dataset
    labels.eagleShuffledByCol = Shuffle(labels.eagleman);
    %     get n matches for each subject
    matches = labels.eagleShuffledByCol == labels.magnet;
    nummatchesinsim = sum(matches,2);
    %     coount number of matches
    temp = hist(nummatchesinsim,0:26);
    %     add to simulation
    sim.colshuffled = sim.colshuffled+temp;
    
    
    
    %     row shuffled data
    % make a simulated dataset
    labels.eagleShuffledByRow = Shuffle(labels.eagleman')';
    %     get n matches for each subject
    matches = labels.eagleShuffledByRow == labels.magnet;
    nummatchesinsim = sum(matches,2);
    %     coount number of matches
    temp = hist(nummatchesinsim,0:26);
    %     add to simulation
    sim.rowshuffled = sim.rowshuffled+temp;
    
    
    
    %     rich distribution
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
    %     get n matches for each subject
    matches = labels.rich == labels.magnet;
    nummatchesinsim = sum(matches,2);
    %     coount number of matches
    temp = hist(nummatchesinsim,0:26);
    %     add to simulation
    sim.rich = sim.rich+temp;
    
    
    % %     uniform distribution across color labels
    % generalte random labels
    labels.uniform = randi(11,size(labels.eagleman))-1;
     %     get n matches for each subject
    matches = labels.uniform == labels.magnet;
    nummatchesinsim = sum(matches,2);
    %     coount number of matches
    temp = hist(nummatchesinsim,0:26);
    %     add to simulation
    sim.uniform = sim.rich+temp;
    
    
    % %    random RGB values
    % % make a set of random rgb values
    rgb.random = rand(size(rgb.eagle));
    %     convert to labels
    labels.random = fpRGB2ColorsJW(rgb.random);
    %     get n matches for each subject
    matches = labels.random == labels.magnet;
    nummatchesinsim = sum(matches,2);
    %     coount number of matches
    temp = hist(nummatchesinsim,0:26);
    %     add to simulation
    sim.randomRGB = sim.rich+temp;
    
    
end


% turn into probability distributions
p.colshuffled = sim.colshuffled/sum(sim.colshuffled);
p.rowshuffled = sim.rowshuffled/sum(sim.rowshuffled);
p.rich = sim.rich/sum(sim.rich);
p.uniform = sim.uniform/sum(sim.uniform);
p.randomRGB = sim.randomRGB/sum(sim.randomRGB);





magnethist = hist(sum(nummatches.eagle,2),0:26);




format short g;

% make a table
disp([(0:26)' fliplr(cumsum(p.rowshuffled(end:-1:1)))' ....
    fliplr(cumsum(p.colshuffled(end:-1:1)))',...
    fliplr(cumsum(p.rich(end:-1:1)))',...
    fliplr(cumsum(p.uniform(end:-1:1)))',...
    fliplr(cumsum(p.randomRGB(end:-1:1)))',...
    fliplr(cumsum(magnethist(end:-1:1)))'])


% figure that shows the subjects sorted by similarity to magnet set with
% the points at which there are n matches labeled

% get index

[y ranking] = sort(nummatches.eagle,'descend');



% 
% % if you wanted to compare it to the earlier figure, you would shuffle the
% % data instead of the magnet set and then sort.  the magnet set is just one
% % throw instead of n?  would be good to use rich, but it is only simulated
% % labels, not really good to simulate rgb values, though we could
% 
% % All the shuffled matches
% % figure('name', ['top ' num2str(length(find(syntype==1))) ' shuffled matches to high fq template'], 'Color', [1 1 1]);
% figure
% 
% [y ranking] = sort(nummatches.eagle,'descend');
% 
% % make a graphical legend 10% as tall as the result
% theLegend = rgb.magnets(1:round(length(find(syntype==1))/10),:,:);
% theResult = rgb.eagle(ranking(100:300),:,:);
% 
% theStack = [theResult; theLegend];
% imagesc(permute(theStack, [1 3 2]));
% 
% set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18);








% lets just count a subset of the column shuffled data






