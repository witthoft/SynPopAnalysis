% script finds various letter by letter cooccurences.  assumes you have
% %
% labels =
%
%               eagleman: [6588x26 double]
%     eagleShuffledByRow: [6588x26 double]
%     eagleShuffledByCol: [6588x26 double]
%                 magnet: [6588x26 double]
%                   rich: [6588x26 double]
%                     fq: [6588x26 double]
%                     random
%                     uniform



eagle2magnetMatches = labels.eagleman == labels.magnet;
% this is a 6588x26 matrix with a 1 everywhere there is a match to the
% magnet set.


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %  are the columns independent of one another?
%   
% suppose that you are exposed to the magnet set.   if you have learned one
% letter, does that increase the likelihood that you learned others?  or is
% it that given the magnet set, the probability that you learn a particular
% letter has some value, and it is independent of wether or not you learn
% any of the other letters.
% we can assess this by testing the independence of the columns.  if two
% events A and B are independent, then the P(A&B)=P(A)*P(B)
makeIndPlots(eagle2magnetMatches,1:26)
% figurename
set(gcf,'Name', 'independence of matches to magnet set, whole data set');







% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % do the same analysis but just on the subset that are magnet
% synesthetes.  want to know if the columns are independent just inside
% that group

% just get the magnet syns
magnetSynsmatches=eagle2magnetMatches(find(syntype==2),:);

makeIndPlots(magnetSynsmatches,1:26);
% figurename
set(gcf,'Name', 'independence of matches to magnet set, just magnet syns');





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % Now do the same thing for frequency



eagle2fqMatches = labels.eagleman == labels.fq;
% this is a 6588x26 matrix with a 1 everywhere there is a match to the
% high frequency template.
% what we want to know is how often is there a one in a given pair of
% columns for a given subject
% you could sum the columns, find the twos, and divide the number by 6588

makeIndPlots(eagle2fqMatches,1:26);
set(gcf,'Name', 'independence of matches to modal choices, all subjects');

% could just look at letters with strong tendencies
makeIndPlots(eagle2fqMatches,[1 2 3 8 9 13 14 15 18 23 24 25 26]);
set(gcf,'Name', 'independence of matches to just letters with strong modal choices, all subjects');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % the column shuffled data is actually independent.  what doees
% that look like
colshuf2fqmatches = labels.eagleShuffledByCol == labels.fq;

makeIndPlots(colshuf2fqmatches,1:26);
set(gcf,'Name', 'independence of matches column shuffled data, all subjects');

% % % % % % % just the subjects with n or more matches

fqSynsmatches=eagle2fqMatches(find(syntype==1),:);
makeIndPlots(fqSynsmatches,1:26);
set(gcf,'Name', 'independence of matches to modal choices, subjects with many matches to modal choice');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % what about the subjects not assigned to a group?
nogroupsynsmatches = eagle2fqMatches(find(syntype==0),:);
makeIndPlots(nogroupsynsmatches,1:26);
set(gcf,'Name', 'independence of matches to modal choices, subjects few matches');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % truly random data sets should also be independent with respect to the
% most frequent or magnet matches
% % 

uniform2fqmatches = labels.uniform == labels.fq;

makeIndPlots(uniform2fqmatches,1:26);
set(gcf,'Name', 'independence of matches to modal choices, uniform distribution, all subs');
% or to the magnet matches
uniform2magmatches = labels.uniform == labels.magnet;

makeIndPlots(uniform2magmatches,1:26);
set(gcf,'Name', 'independence of matches to magnets, uniform distribution, all subs');


% % % % % % % % % % % % % % % % 
randrgb2fqmatches = labels.random == labels.fq;

makeIndPlots(randrgb2fqmatches,1:26);
set(gcf,'Name', 'independence of matches to modal, rand rgb, all subs');

% or to the magnet matches
randrgb2magmatches = labels.random == labels.magnet;

makeIndPlots(randrgb2magmatches,1:26);
set(gcf,'Name', 'independence of matches to magnets, rand rgb, all subs');




% mixture data set
mixture = [magnetSynsmatches; randrgb2magmatches(1:5000,:)];

makeIndPlots(mixture,1:26);
set(gcf,'Name', 'independence of matches to modal, rmagnets mixed with rand rgb');


% try using the subjects with few matches to magnet set
magmatchesindex = sum(eagle2magnetMatches,2);
[y indx] = sort(magmatchesindex);
dohavemagnets = eagle2magnetMatches(indx(1:400),:);
% get worst subjects
indx = flipud(indx);
donthavemagnets = eagle2magnetMatches(indx(1:5000),:);

makeIndPlots([dohavemagnets;donthavemagnets],1:26);


% let's plot each group

% magnets
% matches
figure('name', 'data sorted by similarity to magnet set', 'Color', [1 1 1]);

% make a graphical legend 10% as tall as the result
theLegend = rgb.magnets(1:length(rgb.magnets)/10,:,:);
% get sorted index
[magnetsim Indx] = sort(nummatches.eagle,'descend');



theResult = rgb.eagle(Indx,:,:);

theStack = [theResult; theLegend];
imagesc(permute(theStack, [1 3 2]))

set(gca, 'XTick', 1:26, 'XTickLabel', letters, 'YTick', [1 n], 'FontSize', 18)
title('data sorted by similarity to magnet set')



% let's assume that the modal choices in the population and the colors of
% the magnet set represent competing influences.  we might therefore expect
% to find that the stronger the influence of a modal population choice, the
% less likely the magnet syns are to match the magnet set.  which is what
% we do find

% let's plot them against each other, but separately for shared and not
% shared letters
shared = [1 3 8 13 14 23];
notshared = setdiff(1:26,shared);
magnetSynsmatches=eagle2magnetMatches(find(syntype==2),:);
genpopmatches = eagle2fqMatches(find(syntype~=2),:);

mbaserate = sum(magnetSynsmatches)/length(magnetSynsmatches);
fbaserate = sum(genpopmatches)/length(genpopmatches);
% compute a correlation

figure('Color',[1 1 1]);

[rhoNS pvalNS] = corr(mbaserate(notshared)',fbaserate(notshared)')
subplot(1,2,1)
text(fbaserate(notshared),mbaserate(notshared),letters(notshared));
box off;
xlabel('proportion letter matches color from common template whole pop');
ylabel('proportion letter matches color from magnet set magnet syns');
set(gca,'YLim',[0 .8],'XLim',[0 .8]);
axis square;
title('letters that differ between magnet and common');
% 
% [rhoS pvalS] = corr(mbaserate(shared)',fbaserate(shared)')
% subplot(1,3,3)
% text(fbaserate(shared),mbaserate(shared),letters(shared));
% box off;
% xlabel('proportion letter matches color from common template whole pop');
% ylabel('proportion letter matches color from magnet set magnet syns');
% title('shared letters between magnet and common');
% axis square;
% 
% 

% we might also expect to find that the probability of the modal match in
% the general population is correlated with the probability of a modal
% match in the magnet population
magnetSynsmatchestofq = eagle2fqMatches(find(syntype==2),:);
fbaserateinm = sum(magnetSynsmatchestofq)/length(magnetSynsmatchestofq);
[rho pval] = corr(fbaserate(notshared)',fbaserateinm(notshared)')
subplot(1,2,2)
text(fbaserate(notshared),fbaserateinm(notshared),letters(notshared));
box off;
xlabel('proportion letter matches color from common template whole pop');
ylabel('proportion letter matches color from magnet set magnet syns');
set(gca,'YLim',[0 .6],'XLim',[0 .6]);
axis square;
title('letters that differ between magnet and common');
% 
% plot2svg('competingcontingencies.svg');


% 
% 
% 
% 
% % here is a situation which could hurt our argument.  suppose you really do
% % have two types of synesthetes.  the first has been exposed to some kind
% % of template and they learn each letter with some probability, and those
% % probabilities are independent.   now suppose that the probability that
% % those were never exposed to the template learn those same pairing is low
% % but also independent.   what will the distribution look like then?  I
% % think this is the answer to the weirdness with the magnet set.  since the
% % columns are relatively independent for just the magnet set syns, and
% % likely random across the rest of the population (just with a different
% % frequency) the whole thing looks independent.  given that argument, it
% % might not be possible to say that all the synesthetes have a single
% % probability of learning each letter color combination, since it might be
% % a mixture of probabilities as long as they are always independent.  let's
% % try this by simulation.
% 
% 
% % want to generate a set of columns independently of one another filled
% % with some proportion of 1s.   Then generate a second set of columns with
% % a much smaller proportion of 1s, but which are also independent.
% % presumably all the columns will appear to be independent.
% 
% % strong group
% 
% % make the matrix.  lets say 500 subs with 10 events each
% 
% strong = rand(500,10);
% 
% % set some column probabilities
% 
% colprobs = [.7 .4 .5 .8 .3 .6 .7 .8 .7 .5];
% 
% % make our matrix binary based on the probability
% 
% for i=1:size(strong,2)
%     %    for each column, if bigger than prob make a match
%     m = find(strong(:,i)<=colprobs(i));
%     strong(m,i)=1;
%     strong(setdiff(1:length(strong),m),i)=0;
% end
% 
% 
% 
% 
% % rest of population
% weak = rand(6000,10);
% % set col probabilities
% wcolprobs = .09*ones(1,10);
% 
% % binarize matrix
% for i=1:size(weak,2)
%     %    for each column, if bigger than prob make a match
%     m = find(weak(:,i)<=wcolprobs(i));
%     weak(m,i)=1;
%     weak(setdiff(1:length(weak),m),i)=0;
% end
% 
% % combine the matrices
% simdata = [strong; weak];
% 
% 
% % now test for independence between the columns across the whole population
% 
% jointsim=[];
% % P(A&B)
% for i=1:10
%     for j=1:10
%         jointsim(i,j)=length(find(simdata(:,i)+...
%             simdata(:,j)==2))/length(simdata);
%     end
% end
% 
% 
% % P(A)*P(B)
% % baserate is the diagonal of the joing matrix
% 
% baser = repmat(diag(jointsim)',10,1);
% baser = baser.*baser';
% 
% 
% % P(A&B)-P(A)*P(B)
% % get upper triangular parts
% jointsimV=[];
% baserV=[];
% for i=1:length(jointsim)
%     jointsimV = [jointsimV jointsim(i,find(1:length(jointsim)>i))];
%     baserV = [baserV baser(i,find(1:length(jointsim)>i))];
% end
% 
% 
% % then let's make the figures:
% figure('Name', 'Inependence Analysis ind colums with mixture of rates','Color',[1 1 1]);
% % P(A&B)
% subplot(2,2,1);
% imagesc(jointsim);
% h=colorbar;
% title('P(A&B)');
% set(gca,'CLim',[0 1],'XTick',[1:26],'YTick',[1:26],'XTickLabel',letters,'YTickLabel',letters);
% % axis equal;
% box off;
% 
% % P(A)*P(B)
% subplot(2,2,2);
% imagesc(baser);
% h=colorbar;
% title('P(A)*P(B)');
% set(gca,'CLim',[0 1],'XTick',[1:26],'YTick',[1:26],'XTickLabel',letters,'YTickLabel',letters);
% % axis equal;
% box off;
% 
% 
% % histogram of P(A&B) - P(A)*P(B)
% subplot(2,2,3);
% hist(jointsimV-baserV,[-.2:.025:.2]);
% box off;
% xlabel('P(A&B)-P(A)*P(B)');
% ylabel('number of letters');
% 
% 
% % scatterplot P(A&B) vs P(A)*P(B)
% subplot(2,2,4);
% 
% % x =y
% plot(0:.05:.9,0:.05:.9,'k--');
% axis equal;
% hold on;
% % since the joint and product matrices are symmetric only need the first
% % half
% scatter(baserV,jointsimV,'ro');
% box off;
% set(gca,'XLim',[0 .35],'YLim',[0,.35]);
% ylabel('P(A&B)');
% xlabel('P(A)*P(B)');
% axis equal;
% 
% % now let's generate the matching analysis.
% 
% % need some row shuffled data
% simdatarshuffle = simdata;
% 
% % you can't just shuffle the rows because the number of matches for each
% % subject will remain the same.  might require generating an actual
% % distribution across alternatives.
% 
% for ii = 1:length(simdata)
%   
%    simdatarshuffle(ii,:) = simdata(ii,Shuffle(1:10));
% end
% 
% 
% figure('name','matching analysis on independent sim with fq mixtures','Color',[1 1 1]);
% plot(0:10,hist(sum(simdata,2),0:10),'k--',0:10,hist(sum(simdatarshuffle,2),0:10),'ro');
% box off;
% 
% 
% 
% 
% % thinking about the whole data set and the independence of random
% % variables.  suppose we wanted to search for contingencies in the data.
% 
% % let's say for every pair of letters a and b,we have 11 possible colors
% % with labels 0-10.  you can make a matrix where each row is determined by
% % the color of letter A, and the columns are the probability of the jth
% % color of B.  If the two are independent, these should all look the same.
% 
% % you would have to do this for all 324 possible pairs of columns
% 
% % these could be summarized as correlation plots between the distributions
% 
% 
% 
% %  this is still simpler than considering all possible templates since
% %  there are n different samples where you have 26 ordered elements each
% %  with 11 alternatives
% 
% 










