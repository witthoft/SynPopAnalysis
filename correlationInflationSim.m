% binning the data as is done in the watson and brang synesthesia papers
% inflates the correlations, often by a factor of 2.  at the moment
% I think this is
% because the analysis moves away from a random effects analysis and
% towards a fixed effects analysis without compensating for that.


% in the brang and watson papers they have some predictor x (some similarity
% rating between a pair of letters
% and a response measure y (average distance in color space between a pair of letters)

% for the regression they bin the letters into 65 5 pair groups by
% the similarity between the letters in the color space


% I have to think about whether or not there is a difference between
% binning in x and binning in y


% so let's simulate

% stdevs
s1 = 1;
s2=s1;
% means
m1=0;
m2=m1;
% correlation
p=[0 .2 .4 .4 .5 .6  .8 1];

% colors for pvals
clrs  = [ 1 0 0;
    1 .5 0;
    .5 .5 0;
    0 1 0;
    0 1 .5;
    0 1 1;
    0 .5 1;
    0 0 1;
    .5 0 1;
    1 0 1];



% n
n=325;
% bin
binsize =5;


  % open figure window
    figure('Name', 'correlation inflation for binning dependent variable', 'Color', [1 1 1])
% for each value of p

for j=1:length(p)
    % number of iterations
    numboots = 1000;
    
    % results
    rval = zeros(numboots,1);
    pval = rval;
    bin_rval = pval;
    bin_pval = pval;
    
    
  
    
    for b=1:numboots
        
        % generate random vectors
        u = randn(1,n);
        v = randn(1,n);
        
        %  s1 and s2 are the standard dev, ms are means, and p is correlation
        x = s1*u+m1;
        y = s2*(p(j)*u+sqrt(1-p(j)^2)*v)+m2;

%         y = s2*(p(j)*u+(1-p(j)^2)^(1/8)*v)+m2;
%  compute correlation
        [rval(b) pval(b)]=corr(x',y','type','Spearman');
        
        %  bin data by sorting y
        
        [sortedy, sortedindx] = sort(y);
        
        % average our responses
        % reshape for averaging
        sortedy=reshape(sortedy,binsize,n/binsize);
        % average
        ave_sortedy=mean(sortedy);
        
        % sort and bin our x values for this boot
        sortedx = reshape(x(sortedindx),binsize,n/binsize);
        ave_sortedx=mean(sortedx);
        
        % compute correlation on binned values
        [bin_rval(b), bin_pval(b)]=corr(ave_sortedx',ave_sortedy','type','Spearman');
        
    end
    
    %add to plots
    subplot(1,2,1);
    hold on;
    scatter(rval.^2,bin_rval.^2,10,clrs(j,:));
    
    
    subplot(1,2,2);
    hold on;
    scatter(pval,bin_pval,10,clrs(j,:));
    
    
end

subplot(1,2,1);
hold on;
xlabel('unbinned r^2');
ylabel('binned r^2');
box off;
set(gca,'Xlim',[0 1],'YLim',[0 1]);
axis square;
plot(0:01:1,0:01:1,'k-');



subplot(1,2,2);
hold on;
xlabel('unbinned p');
ylabel('binned p');
box off;
set(gca,'Xlim',[0 1],'YLim',[0 1]);
axis square;
% axis equal;
plot(0:0.05:1,0:0.05:1,'k-')

% in general we find that the r^2 is always inflated while the p-value does
% not change.  is this an error?  don't have to look at it very long to
% figure out that you could analytically determine the amount of expected
% inflation of variance explained for smoothing given the parameters.  

% depends on how you interpret it I suppose, if you think the variance you
% are explaining applies to particular letter pairs, then no.  and that
% might be the thing you care about.

% its also not clear to me that the p-value makes sense....
%  you are reducing your n which is why the
% p-value goes down even though the correlations go up.
% so in terms of significance, your significant effects don't really
% change.  but the assumption is that all of your datapoints are
% independent samples and this isn't true any more?  if not why not?



% for kicks do as differences.  would have been smarter to save the
% variables than copy and past

  % open figure window
    figure('Name', 'correlation inflation for binning dependent variable', 'Color', [1 1 1])
% for each value of p

for j=1:length(p)
    % number of iterations
    numboots = 1000;
    
    % results
    rval = zeros(numboots,1);
    pval = rval;
    bin_rval = pval;
    bin_pval = pval;
    
    
  
    
    for b=1:numboots
        
        % generate random vectors
        u = randn(1,n);
        v = randn(1,n);
        
        %  s1 and s2 are the standard dev, ms are means, and p is correlation
        x = s1*u+m1;
        y = s2*(p(j)*u+sqrt(1-p(j)^2)*v)+m2;
        %  compute correlation
        [rval(b) pval(b)]=corr(x',y','type','Spearman');
        
        %  bin data by sorting y
        
        [sortedy, sortedindx] = sort(y);
        
        % average our responses
        % reshape for averaging
        sortedy=reshape(sortedy,binsize,n/binsize);
        % average
        ave_sortedy=mean(sortedy);
        
        % sort and bin our x values for this boot
        sortedx = reshape(x(sortedindx),binsize,n/binsize);
        ave_sortedx=mean(sortedx);
        
        % compute correlation on binned values
        [bin_rval(b), bin_pval(b)]=corr(ave_sortedx',ave_sortedy','type','Spearman');
        
    end
    
    %add to plots
    subplot(2,1,1);
    hold on;
    plot(rval.^2,bin_rval.^2-rval.^2,'o','MarkerEdgeColor',clrs(j,:),...
        'MarkerFaceColor',clrs(j,:),'MarkerSize',1);
    
    
    subplot(2,1,2);
    hold on;
    plot(pval,bin_pval-rval,'o','MarkerEdgeColor',clrs(j,:),...
        'MarkerFaceColor',clrs(j,:),'MarkerSize',1);
    
    
end

subplot(2,1,1);
hold on;
xlabel('unbinned r^2');
ylabel('ubinned - binned r^2');
box off;
set(gca,'Xlim',[0 1],'YLim',[-1 1]);




subplot(2,1,2);
hold on;
xlabel('unbinned p');
ylabel('unbinned - binned p');
box off;
set(gca,'Xlim',[0 .1],'XScale','log','YLim',[-1 1]);



% even given all this, it is still a kind of fixed effects analysis because
% you are averaging out the subject variation.  if your subjects are
% allowed to eat up some variance...


