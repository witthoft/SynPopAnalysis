
% thinking about the whole data set and the independence of random
% variables.  suppose we wanted to search for contingencies in the data.

% let's say for every pair of letters a and b,we have 11 possible colors
% with labels 0-10.  you can make a matrix where each row is determined by
% the color of letter A, and the columns are the probability of the jth
% color of B.  If the two are independent, these should all look the same.

% you would have to do this for all 324 possible pairs of columns


% set the labels
syntype = zeros(length(dbLabeled),1);


% num of matches needed to be in group
magnetthreshold = 10;
fqthreshold = 10;

% code for labels variable is in s_fpEaglemanMatches but in this script
% should be in workspacelabels.
magmatches = labels.eagleman == labels.magnet;
fqmatches = labels.eagleman == labels.fq;
% let's set the fq first
syntype(find(sum(fqmatches,2)>=fqthreshold))=1;
% then magnet.
syntype(find(sum(magmatches,2)>=magnetthreshold))=2;
% maybe want to find those that are in both groups?





% which labels to work on
% % database
% whichlabels = labels.eagleman; %database
% savedir = 'indRandVarPlots/database/';
% saveprefix = 'eagleman';

% % column shuffled
% whichlabels = labels.eagleShuffledByCol ;%column shuffled labels
% savedir = 'indRandVarPlots/colshuffled/';
% saveprefix = 'colshuffled.';

% % % uniform
% whichlabels = labels.uniform;
% savedir = 'indRandVarPlots/uniform/';
% saveprefix = 'uniform.';


% just the magnet synesthetes (more than 10 matches)
whichlabels = labels.eagleman(find(syntype==2),:);
savedir = 'indRandVarPlots/magnetsyns/';
saveprefix = 'magnetsyns.';








% 
% % all but the magnet synesthetes
% whichlabels = labels.eagleman(find(syntype~=2),:);
% savedir = 'indRandVarPlots/allbutmagnets/';
% saveprefix = 'notmagnetsyns.';



% these could be summarized as correlation plots between the distributions

% check for directory
if ~isdir(savedir)
    mkdir(savedir)
end


% data is labels.eagleman
colornames = {...
    'black' ...
    'white' ...
    'red' ...
    'green' ...
    'yellow' ...
    'blue' ...
    'brown' ...
    'purple' ...
    'pink' ...
    'orange' ...
    'grey'...
    'none'
    };
%
% for each letter
for lA=1:26
%     make a variable to store the probability distributions
    BcA =zeros(12);
    %     pick another letter (only need upper triangle
    for lB=1:26
        %     for each color category
        for cA=1:12
            %     get a histogram of the values of lB
            %     get index to entries where letter A has color a
            lAcA = find(whichlabels(:,lA)==cA-1);
            %          then get all of the values of B for those subjects
            BcA(cA,:) = hist(whichlabels(lAcA,lB),0:11);
%             turn this into a probability distribution?
            BcA(cA,:) = BcA(cA,:)/sum(BcA(cA,:));
        end
%         make a figure with the histogram comparing the letters
        figure('Name',['joint distribution of ' letters(lA) 'and ' letters(lB)],...
            'Color',[1 1 1]);
        imagesc(BcA);
        box off;
        set(gca,'XTick',1:12,'XTickLabel',colornames,'YTick',1:12,'YTickLabel',colornames);
        xlabel(['probability ' letters(lB) ' is ']);
        ylabel(['if ' letters(lA) ' is']);
        set(gca,'CLim',[0 .5]);
        colorbar;
%         save the figure
        saveas(gcf,[savedir saveprefix letters(lA) 'vs' letters(lB) '.png'],'png');
        plot2svg([savedir saveprefix letters(lA) 'vs' letters(lB) '.svg'],gcf);
        close(gcf);

    end
end




% make a large single figure with all the letters vs all the letters so
%make indices so we know where to place each probability distribution we
%solve for
xindvector = ((0:26)*12)+1;
% matrix to hold our solutions
BcA = zeros(12*26);
% for each letter
for lA=1:26
    %     make a variable to store the probability distributions
    
    %     pick another letter (only need upper triangle
    for lB=1:26
        %     for each color category
        for cA=1:12
            %     get a histogram of the values of lB
            %     get index to entries where letter A has color a
            lAcA = find(whichlabels(:,lA)==cA-1);
            %          then get all of the values of B for those subjects
            temp = hist(whichlabels(lAcA,lB),0:11);
            %             turn this into a probability distribution?
            temp = temp/sum(temp);
            %             inset into our matrix
            BcA(xindvector(lA)+cA-1,xindvector(lB):xindvector(lB)+11)=temp;
        end
        
    end
end

% make a figure with the histogram comparing the letters
figure('Name',['joint distribution of ' letters(lA) 'and ' letters(lB)],...
    'Color',[1 1 1]);
imagesc(BcA);
box off;
set(gca,'XTick',6:12:306,'XTickLabel',letters,'YTick',6:12:306,'YTickLabel',letters);
xlabel(['probability ' letters(lB) ' is ']);
ylabel(['if ' letters(lA) ' is']);
set(gca,'CLim',[0 .5]);
colorbar;
%         save the figure
saveas(gcf,[savedir saveprefix 'fullmatrix.jpg'],'jpg');
plot2svg([savedir saveprefix 'fullmatrix.svg'],gcf);

close(gcf);

% number of iterations of model fits
modelits = 50;



% suppose we want to model the conditional probability structure
% we can split our data in half and use one half of the data to model the
% other half

% models

% color matching is random
% here our model is uniform color selection for each letter and the letters
% are independent of one another

% test data is just the dataset
    testdata =whichlabels;

    CPMmodelAll=zeros(312);

for i=1:modelits
    
    modeldata = randi(12,size(testdata))-1;
    %      get the conditional probability matrix for both
    CPMmodel = getCondPMatrix(modeldata);
    CPMmodelV = reshape(CPMmodel,length(CPMmodel)^2 ,1);
    UniformCPMtest = getCondPMatrix(testdata);
    UniformCPMtestV = reshape(UniformCPMtest,length(UniformCPMtest)^2 ,1);
    %     unpack rows into a vector
    %      compute a correlation
  
    [uniformr(i), uniformp(i)]  = corr(CPMmodelV,UniformCPMtestV,'rows','pairwise') ;
    
    CPMmodelAll = CPMmodelAll+CPMmodel;
    
end

% r = .57 on average

CPMmodelAll = CPMmodelAll/modelits;

% make a figure with the histogram comparing the letters
figure('Name',['data generated using uniform model'],...
    'Color',[1 1 1]);
imagesc(CPMmodelAll);
box off;
set(gca,'XTick',6:12:306,'XTickLabel',letters,'YTick',6:12:306,'YTickLabel',letters);
xlabel(['probability ' letters(lB) ' is ']);
ylabel(['if ' letters(lA) ' is']);
set(gca,'CLim',[0 .5]);
colorbar;













% color matching is biased for each letter, but letters are independent
%  model is to build a sample data set by using the marginal distribution
%  for each letter

% so need to generate subjects by generating data only knowing the base
% rates for each letter-color assignment

% want to get a histogram of the right size for each column.  then turn
% that into a probability distribution

% get counts for labels in each column
[letterlabelcounts X] =hist(whichlabels,0:11);

% convert those to probabilities
letterlabelprobs = letterlabelcounts/length(whichlabels);

% check
sum(letterlabelprobs)


% now build a sample using the empirical marginal distributions



% or you could just generate random indices into the columns.  this is just
% column shuffling the data. 

% 
columnMarginModel = getCondPMatrix(labels.eagleShuffledByCol);

columnMarginModelV = reshape(columnMarginModel,length(columnMarginModel)^2,1);

[cmmr(1), cmmpp(1)]  = corr(columnMarginModelV,CPMtestV,'rows','pairwise') ;

 % make a figure with the histogram comparing the letters
figure('Name',['data generated using colum probabilities'],...
    'Color',[1 1 1]);
imagesc(columnMarginModel);
box off;
set(gca,'XTick',6:12:306,'XTickLabel',letters,'YTick',6:12:306,'YTickLabel',letters);
xlabel(['probability ' letters(lB) ' is ']);
ylabel(['if ' letters(lA) ' is']);
set(gca,'CLim',[0 .5]);
colorbar;





    SplitHalfmodelAll=zeros(312);


% split half
% given any half of the data, how similar is it to the other half
for i=1:modelits
    %     get random order of subjects
    x=Shuffle(1:length(whichlabels));
    modeldata = x(1:length(whichlabels)/2);
    testdata = x(length(whichlabels)/2+1:end);
    %      get the conditional probability matrix for both
    CPMmodel = getCondPMatrix(whichlabels(modeldata,:));
    CPMmodelV = reshape(CPMmodel,length(CPMmodel)^2 ,1);
    CPMtest = getCondPMatrix(whichlabels(testdata,:));
    CPMtestV = reshape(CPMtest,length(CPMtest)^2 ,1);
    %     unpack rows into a vector
    %      compute a correlation
  
    [splithalfr(i), splithalfp(i)]  = corr(CPMmodelV,CPMtestV,'rows','pairwise') ;
    
end
% .93 on average.  very tight distribution of rs

%






%
% % let's code it again for sanity purposes.  it must be right but there are
% % a few confusing things about it.
%
% % the task is given a set of subjects each one of 12 labels for each
% % letter, generate a matrix of the following type,  given a letter and a
% % color, what is the probability for another letter to have each of the
% % colors.
% % matrix to hold our random variable space
% AvsB = zeros(26*12);
% whichrow = 1;
% % for each letter
% for letA=1:26
% %     pick a color category
%     for colorcatA=0:11
% %         find all the entries in A that correspond to that color
%             letAandColA = find(whichlabels(:,letA)==colorcatA);
% %         fill in that whole row in the matrix
%             letterArow = [];
% %         pick a second letter
%         for letB=1:26
% %             find distribution across colors
%             letBdist = hist(whichlabels(letAandColA,letB),0:11);
% %             turn into probability dist
%             letBdist = letBdist/sum(letBdist);
%             letterArow = [letterArow letBdist];
%         end
% %         add new row to matrix
%         AvsB(whichrow,:)=letterArow;
%         whichrow=whichrow+1;
%     end
% end
%
% % make a figure with the histogram comparing the letters
% figure('Name',['conditional distribution of all letters'],...
%     'Color',[1 1 1]);
% imagesc(AvsB);
% box off;
% set(gca,'XTick',6:12:306,'XTickLabel',letters,'YTick',6:12:306,'YTickLabel',letters);
% xlabel(['given ']);
% ylabel(['prob across']);
% set(gca,'CLim',[0 .5]);
% colorbar;
%
% % save the full probability matrix into a .mat variable
% save(['indRandVarPlots/' saveprefix 'fullmatrix.mat'],'AvsB');
%
%

%
%
%
%
% % what about the joint probability matrix?
%
%
% % let's code it again for sanity purposes.  it must be right but there are
% % a few confusing things about it.
%
% % the task is given a set of subjects each one of 12 labels for each
% % letter, generate a matrix of the following type,  given a letter and a
% % color, what is the probability for another letter to have each of the
% % colors.
% % matrix to hold our random variable space
% AvsB = zeros(26*12);
% whichrow = 1;
% % for each letter
% for letA=1:26
% %     pick a color category
%     for colorcatA=0:11
% %         find all the entries in A that correspond to that color
%             letAandColA = find(whichlabels(:,letA)==colorcatA);
% %         fill in that whole row in the matrix
%             letterArow = [];
% %         pick a second letter
%         for letB=1:26
% %             find distribution across colors
%             letBdist = hist(whichlabels(letAandColA,letB),0:11);
% %             turn into probability dist
%             letBdist = letBdist/sum(letBdist);
% %             multiply by the probability that letA is colorcatA
%             letBdist = letBdist*(length(letAandColA)/length(whichlabels));
%             letterArow = [letterArow letBdist];
%         end
% %         add new row to matrix
%         AvsB(whichrow,:)=letterArow;
%         whichrow=whichrow+1;
%     end
% end
%
%
%
% % make a figure with the histogram comparing the letters
% figure('Name',['joint distribution of all letters'],...
%     'Color',[1 1 1]);
% imagesc(AvsB);
% box off;
% set(gca,'XTick',6:12:306,'XTickLabel',letters,'YTick',6:12:306,'YTickLabel',letters);
% xlabel(['given ']);
% ylabel(['prob across']);
% set(gca,'CLim',[0 .1]);
% colorbar;
%
% % save the full probability matrix into a .mat variable
% save(['indRandVarPlots/' saveprefix 'fullmatrix.mat'],'AvsB');
%
%
%
%










% to double check , would be niced to be able to do a column

%  this is still simpler than considering all possible templates since
%  there are n different samples where you have 26 ordered elements each
%  with 11 alternatives


% do PCA on the full matrix
% what are we asking?
% we are choosing a point wheich is a letter+color assignment.  there are
% 28*11 of these points.  each point lies somewhere in a 28*11 dimensional
% space, where the location on each dimension is the p that another letter
% had a particular color letter assignment.  so points are things like, a is red.  and
% dimensions are also things like a is red, or b is blue.  we can imagine
% trying to reduce the dimensionality of the data by rotating the axes.
% the values on each dimension only go from 0 to 1... does it matter that
% the data has been transformed from counts to probabilities?  a weird
% feature is that any point is at one on the dimension that is the same as
% its label (i.e. a is red is at 1 on the dimension a is red).
%
% %do pca
%     [pca.pcs pca.score pca.latent pca.tsquare] = princomp(BcA);
%
%     %find amount of variance explained by adding additional pcs
%     pca.varExp = cumsum(pca.latent)/sum(pca.latent);
%
%     %get mean of the original condition space
%     pca.mean = mean(pca.voxels(:,conds));
%
%
%







