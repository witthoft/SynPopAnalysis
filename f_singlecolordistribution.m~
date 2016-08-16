function   colorinfo = f_singlecolordistribution(rgb_data, group, colorcat)

%  colorinfo = singlecolordistribution(rgb_data, group)
% rgb_data is subset of u_rgb
% group is a string used for figure labels designating subset of
% synesthetes used
% output
% take a set of rgb matching data and count the number of letters not
% matched to a color category in various ways (1-12)
% total number of possible matches
% number not matched
% number not matched by letter
% number not matched by subject
% 

names = {...
    'black' ... %0
    'white' ...%1
    'red' ...%2
    'green' ...%3
    'yellow' ...%4
    'blue' ...%5
    'brown' ...%6
    'purple' ...%7
    'pink' ...%8
    'orange' ...%9
    'grey'...%10
    'none',...%nan
    };


clrs = [ 0 0 0;
    1 1 1;
    1 0 0;
    0 1 0;
    1 1 0;
    0 0 1;
    .6 .3 0;
    .4 0 .8;
    1 .2 1;
    1 .5 0;
    .5 .5 .5;
    0 0 0;
    ];


% this is handy for labeling figures
letters = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'};





% frequency of letters as first letter of word in English text





% need to load up dbNumbered
[dbNumbered dbLabeled] = fpRGB2ColorsJW(permute(rgb_data, [1,3,2]));
% [dbNumbered dbLabeled] = fpRGB2ColorsJW(rgb_data);


% index to subjects with colorcategory
subswcolor = [];
% how many of color category for each subject
numcolors = [];
% for each subject
for i=1:length(dbNumbered)
%     check for a colorcat
    y=length(find(dbNumbered(i,:)==colorcat));
    if y>0
%     fprintf('subject %g has %g nans \n',i,y);
    subswcolor = [subswcolor i];
    numcolors = [numcolors y];
    else
        numcolors = [numcolors 0];
    end
end
    


% make a directory to save the figures if one doesn't exist

colorsavdir = 'singlecolordistribution';

if ~exist(colorsavdir,'dir')
    mkdir(colorsavdir);
end


% number of histogram of subjects with n NaNs

figure('Name',['histogram of subjects with n letters with color set to ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));
hist(numcolors,[0:26]);
xlabel('number of letters with match');
% so many are 0 that we need a log scale on y
hold on;
ylabel('number of subjects');
box off;
set(gca,'XLim',[-.5 26]);
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.histofmatchesAllsubs.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.histofmatchesAllsubs.svg'],gcf);

close(gcf);



% do this as percentages
h=hist(numcolors,[0:26]);

pctcolorcat = h/sum(h);




figure('Name',['percent of subjects with n letters matched to ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));
bar(pctcolorcat);

xlabel(['number of letters matched to ' names{colorcat+1}]);
ylabel('percent of subjects');
box off;
set(gca,'XLim',[0 26]);

saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.normhistofnomatchesAllsubs.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.normhistofnomatchesAllsubs.svg'],gcf);
close(gcf);





% how is this behavior distributed across letters?
figure('Name','distribution of nans across letters','Color',[1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,2,1);
% figure where nans are black and color matches are white
% get one rgb column of data for subs x letters
colorcatmatrix = dbNumbered;


% this will behave weirdly for black which is also 0.  what a hack!
if colorcat == 0
    %     set not black to 2
    colorcatmatrix(find(colorcatmatrix~=colorcat))=2;
    colorcatmatrix(find(colorcatmatrix==colorcat))=1;
    colorcatmatrix(find(colorcatmatrix==2))=0;
    
else
    colorcatmatrix(find(colorcatmatrix~=colorcat))=0;
    colorcatmatrix(find(colorcatmatrix==colorcat))=1;
end


imagesc(colorcatmatrix);  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(colorcatmatrix)/length(dbNumbered));
box off;
xlabel('letters');
ylabel('% of times not matched');
set(gca,'XTick',[1:26],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.distofnonmatchesAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.distofnonmatchesAcrossLetters.svg'],gcf);
close(gcf);




% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));

scatter(1:26,sum(colorcatmatrix)/length(dbNumbered),100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(1.25:26.25,sum(colorcatmatrix)/length(dbNumbered),letters,'FontWeight','Bold','FontSize',20);
box off;
xlabel('letters');
ylabel('% of times matched');
set(gca,'XTick',[1:26],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(1:26,sum(colorcatmatrix)/length(dbNumbered),1);
    %   points predicted by fit line
    yfit=polyval(p,1:26);
    %   get the residuals
    yresid=(sum(colorcatmatrix)/length(dbNumbered))-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(sum(colorcatmatrix)/length(dbNumbered));
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((1:26)',(sum(colorcatmatrix)/length(dbNumbered))','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((1:26)',(sum(colorcatmatrix)/length(dbNumbered))','type','Kendall');

% put the regression line on the plot
xpts =1:.25:26;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.scatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.scatterAcrossLetters.svg'],gcf);
close(gcf);
% 
% 
% 







% 
% need to add lewand frequency here if this is going to work. is still in
% main script
% 
% % hmmmm is number of nans predicted by the letter frequency?
% 
% lewand frequency
% 5 20  1 15  9 14 19  8 18  4 12 3 21  13 23 6 7 25  16  2 22   11 10 24 17 26
% e t   a o   i n  s   h r   d l  c u   m  w  f g y p   b v    k  j  x  q  z
fqorder = [5 20  1 15  9 14 19  8 18  4 12 3 21  13 23 6 7 25  16  2 22   11 10 24 17 26];





% how is this behavior distributed across letters?
figure('Name',['distribution of letters matched to ' names{colorcat+1} ' fqsort'],'Color',[1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,2,1);
% figure where nans are black and color matches are white
% get one rgb column of data for subs x letters

imagesc(colorcatmatrix(:,fqorder));  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters(fqorder));
% set(gca,'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(colorcatmatrix(:,fqorder))/length(dbNumbered));
box off;
xlabel('letters');
ylabel(['% of times matched to ' names{colorcat+1}]);
set(gca,'XTick',[1:26],'XTickLabel',letters(fqorder));
% set(gca,'XTickLabel',letters);
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.distofnonmatchesAcrossLettersFQ.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.distofnonmatchesAcrossLettersFQ.svg'],gcf);
close(gcf);




% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));

yfqorder=sum(colorcatmatrix)/length(dbNumbered);
yfqorder=yfqorder(fqorder);

scatter(1:26,yfqorder,100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(1.25:26.25,yfqorder,letters(fqorder),'FontWeight','Bold','FontSize',20);
box off;
xlabel('letters');
ylabel('% of times matched');
set(gca,'XTick',[1:26],'XTickLabel',letters(fqorder));
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(1:26,yfqorder,1);
    %   points predicted by fit line
    yfit=polyval(p,1:26);
    %   get the residuals
    yresid=(sum(yfqorder))-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(yfqorder);
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((1:26)',(yfqorder)','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((1:26)',(yfqorder)','type','Kendall');

% put the regression line on the plot
xpts =1:.25:26;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.firstletterscatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.firstletterscatterAcrossLetters.svg'],gcf);
close(gcf);





% instead of using frequency to set rank use actual frequency

% frequency of letters in English text in alphabetical order

lewandfq = [
0.08167
0.01492
0.02782
0.04253
0.12702
0.02228
0.02015
0.06094
0.06966
0.00153
0.00772
0.04025
0.02406
0.06749
0.07507
0.01929
0.00095
0.05987
0.06327
0.09056
0.02758
0.00978
0.0236
0.0015
0.01974
0.00074];




% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));


scatter(lewandfq,sum(colorcatmatrix)/length(dbNumbered),100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(lewandfq+.0025,sum(colorcatmatrix)/length(dbNumbered),letters,'FontWeight','Bold','FontSize',20);
box off;
xlabel('lewand frequency in English text');
ylabel('% of times not matched');
% set(gca,'XTick',[lewandfq],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(lewandfq',sum(colorcatmatrix)/length(dbNumbered),1);
    %   points predicted by fit line
    yfit=polyval(p,lewandfq);
    %   get the residuals
    yresid=(sum(colorcatmatrix)/length(dbNumbered))'-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(sum(colorcatmatrix)/length(dbNumbered));
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((lewandfq),(sum(colorcatmatrix)/length(dbNumbered))','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((lewandfq),(sum(colorcatmatrix)/length(dbNumbered))','type','Kendall');

% put the regression line on the plot
xpts =0:.005:.14;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.lewandfqletterscatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.lewandfqletterscatterAcrossLetters.svg'],gcf);
close(gcf);





% do it as log lewand fq
lewandfq = log(lewandfq);
% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));


scatter(lewandfq,sum(colorcatmatrix)/length(dbNumbered),100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(lewandfq+.0025,sum(colorcatmatrix)/length(dbNumbered),letters,'FontWeight','Bold','FontSize',20);
box off;
xlabel('log lewand frequency in English text');
ylabel('% of times not matched');
% set(gca,'XTick',[lewandfq],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(lewandfq',sum(colorcatmatrix)/length(dbNumbered),1);
    %   points predicted by fit line
    yfit=polyval(p,lewandfq);
    %   get the residuals
    yresid=(sum(colorcatmatrix)/length(dbNumbered))'-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(sum(colorcatmatrix)/length(dbNumbered));
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((lewandfq),(sum(colorcatmatrix)/length(dbNumbered))','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((lewandfq),(sum(colorcatmatrix)/length(dbNumbered))','type','Kendall');

% put the regression line on the plot
xpts =0:-.25:-8;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.loglewandfqletterscatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.loglewandfqletterscatterAcrossLetters.svg'],gcf);
close(gcf);








% fq of first letter in a word

% 20 1 19 8 23 9 15 2 13 6 3 12 4 16 14 5 7 18 25 21 22 10 11 24 26 24
% t  a s  h w  i o  b m  f c l  d p  n  e g r  y  u  v  j  k  q  z  x

fqfirst = [20 1 19 8 23 9 15 2 13 6 3 12 4 16 14 5 7 18 25 21 22 10 11 17 26 24];

% how is this behavior distributed across letters?
figure('Name','distribution of nans across letters fqsort','Color',[1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,2,1);

% resort x axis
imagesc(colorcatmatrix(:,fqfirst));  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters(fqfirst));
% set(gca,'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(colorcatmatrix(:,fqfirst))/length(rgb_data));
box off;
xlabel('letters');
ylabel('% of times not matched');
set(gca,'XTick',[1:26],'XTickLabel',letters(fqfirst));
% set(gca,'XTickLabel',letters);
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.distofnonmatchesAcrossLettersFQfirst.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.distofnonmatchesAcrossLettersFQfirst.svg'],gcf);
close(gcf);




% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));

yfqfirst=sum(colorcatmatrix)/length(dbNumbered);
yfqfirst=yfqfirst(fqfirst);

scatter(1:26,yfqfirst,100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(1.25:26.25,yfqfirst,letters(fqfirst),'FontWeight','Bold','FontSize',20);
box off;
xlabel('letters');
ylabel('% of times matched');
set(gca,'XTick',[1:26],'XTickLabel',letters(fqfirst));
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(1:26,yfqfirst,1);
    %   points predicted by fit line
    yfit=polyval(p,1:26);
    %   get the residuals
    yresid=(yfqfirst)-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(yfqfirst);
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((1:26)',(yfqfirst)','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((1:26)',(yfqfirst)','type','Kendall');

% put the regression line on the plot
xpts =1:.25:26;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.firstscatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.firstscatterAcrossLetters.svg'],gcf);
close(gcf);







% instead of using frequency of first letter to set rank use actual frequency

% frequency of first letter in words in English text in alphabetical order
fqfirstletter = [
.11602% 	 
.04702% 	
.03511% 	
.02670% 	
.02007% 	
.03779% 	
.01950% 	
.07232% 	
.06286% 	
.00597% 	
.00590% 	
.02705% 	
.04383% 	
.02365% 	
.06264% 	
.02545% 	
.00173% 	
.01653% 	
.07755% 	
.16671% 	
.01487% 	
.00649% 	
.06753% 	
.00017% 	
.01620% 	
.00034];

% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));


scatter(fqfirstletter,sum(colorcatmatrix)/length(dbNumbered),100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(fqfirstletter+.0025,sum(colorcatmatrix)/length(dbNumbered),letters,'FontWeight','Bold','FontSize',20);
box off;
xlabel('frequency first letter in English text');
ylabel('% of times matched');
% set(gca,'XTick',[fqfirstletter],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(fqfirstletter',sum(colorcatmatrix)/length(dbNumbered),1);
    %   points predicted by fit line
    yfit=polyval(p,fqfirstletter);
    %   get the residuals
    yresid=(sum(colorcatmatrix)/length(dbNumbered))'-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(sum(colorcatmatrix)/length(dbNumbered));
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((fqfirstletter),(sum(colorcatmatrix)/length(dbNumbered))','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((fqfirstletter),(sum(colorcatmatrix)/length(dbNumbered))','type','Kendall');

% put the regression line on the plot
xpts =0:.005:.17;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.fqfirstletterletterscatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.fqfirstletterletterscatterAcrossLetters.svg'],gcf);
close(gcf);






% do as log fq of first letter
fqfirstletter=log(fqfirstletter);
% as scatterplot with correlation
figure('Name',['scatterplot of letters x ' names{colorcat+1}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));


scatter(fqfirstletter,sum(colorcatmatrix)/length(dbNumbered),100,'MarkerFaceColor',clrs(colorcat+1,:),...
    'MarkerEdgeColor',[0 0 0]);
hold on
text(fqfirstletter+.0025,sum(colorcatmatrix)/length(dbNumbered),letters,'FontWeight','Bold','FontSize',20);
box off;
xlabel('frequency first letter in English text');
ylabel('% of times matched');
% set(gca,'XTick',[fqfirstletter],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);
% add correlation and regression line

 %     compute a correlation between the measures
    %  linear regression
    p=polyfit(fqfirstletter',sum(colorcatmatrix)/length(dbNumbered),1);
    %   points predicted by fit line
    yfit=polyval(p,fqfirstletter);
    %   get the residuals
    yresid=(sum(colorcatmatrix)/length(dbNumbered))'-yfit;
    %    sumsquaredresidual
    SSresid=sum(yresid.^2);
    %     total sum of squares (num letters-1)*var
    SSTotal=25*var(sum(colorcatmatrix)/length(dbNumbered));
    %     rsquared
    rsq = 1 - SSresid/SSTotal;
    %
    %   compare to matlab output
    [s_rho, pval] = corr((fqfirstletter),(sum(colorcatmatrix)/length(dbNumbered))','type','Spearman');
    
%     let's get a nonparametric rank correlation (tau)
        [s_tau, kpval] = corr((fqfirstletter),(sum(colorcatmatrix)/length(dbNumbered))','type','Kendall');

% put the regression line on the plot
xpts =0:-.25:-8;
plot(xpts,polyval(p,xpts),'k--');

% write down the numbers in the title
title(['r^2 = ' num2str(rsq) ' p = ' num2str(pval)]);


% 
saveas(gcf,[colorsavdir '/' group '.' names{colorcat+1} '.logfqfirstletterletterscatterAcrossLetters.png'],'png');
plot2svg([colorsavdir '/' group '.' names{colorcat+1} '.logfqfirstletterletterscatterAcrossLetters.svg'],gcf);
close(gcf);






% make some tables and return some info

colorinfo.totalmatches = length(dbNumbered)*26;
colorinfo.matchtocategory = sum(numcolors);
colorinfo.catbyletter = sum(colorcatmatrix);
colorinfo.catbysubject = hist(numcolors,0:26);
colorinfo.cat = names{colorcat+1};


% % % % % % % % % % % % % % % % % % 
% % can rerun if needed
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % 


% 
% 
% % how big does sample need to be to estimate the proportion of times a
% % letter is a particular color?   need to do a bootstrap across sample
% % sizes.  let's start with a resolution of 25 and 100 bootstraps for each
% % sample size.  we can put the resulting curve into a series of subplots
% 
% % bootstrap parameters:
% btresolution = 25;
% btsamplesizes = [btresolution:btresolution:length(rgb_data)];
% numboots = 100;
% 
% bootoutputs = nan(numboots,26, length(btsamplesizes));
% 
% % figure
% 
% 
%     %    for each each bootstrap sample size
%     for samplesize=1:length(btsamplesizes)
%         %         do numboots estimates
%         for bt = 1:numboots
% %             get a bootstrap sample
%             bsample = randi(length(rgb_data),[1,btsamplesizes(samplesize)]);
% % count colorcat occurence for bootstrap
%             bootoutputs(bt,:,samplesize)=sum(colorcatmatrix(bsample,:))/btsamplesizes(samplesize);
%             
%         end
%         
%     end
%     
% 
%     
% % want to summarize with median and 95% confidence interval
%     
% ptile = prctile(bootoutputs,[2.5,50,97.5]);
%     
% figure('Name',['Bootstrap confidence interval on probability letter is ' names{colorcat+1}],...
%     'Color','White','Position',get(0,'ScreenSize'));
% 
% for f=1:26    
%     
%     subplot(6,5,f);
%     hold on;
%     plot(btsamplesizes,squeeze(ptile(:,f,:)));
%     ylabel(['p ' letters(f) ' is ' names{colorcat+1}]);
%     xlabel('bootstrap sample size');
%     set(gca,'YLim',[0 .7],'XLim',[0 max(btsamplesizes)]);
% % bootoutput is a 3d matrix bootstrap x letter x sample size
% % plot(btsamplesizes,squeeze(bootoutputs(:,f,:)));
% end
% 
% saveas(gcf,[colorsavdir '/' group '.bootstrapforeachletter.' names{colorcat+1} '.distofnonmatchesAcrossLetters.png'],'png');
% plot2svg([colorsavdir '/' group '.bootstrapforeachletter.' names{colorcat+1} '.distofnonmatchesAcrossLetters.svg'],gcf);
% close(gcf);
% 
%     
% figure('Name',['width of ootstrap confidence interval on probability letter is ' names{colorcat+1}],...
%     'Color','White','Position',get(0,'ScreenSize'));
% 
% for f=1:26    
%     
%     subplot(6,5,f);
%     hold on;
%     plot(btsamplesizes,squeeze(ptile(3,f,:))-squeeze(ptile(1,f,:)));
%     ylabel(['p ' letters(f) ' is ' names{colorcat+1}]);
%     xlabel('bootstrap sample size');
%     set(gca,'YLim',[0 .1],'XLim',[0 max(btsamplesizes)]);
% % bootoutput is a 3d matrix bootstrap x letter x sample size
% % plot(btsamplesizes,squeeze(bootoutputs(:,f,:)));
% end
% 
% 
% saveas(gcf,[colorsavdir '/' group '.bootstrapciwidthforeachletter.' names{colorcat+1} '.distofnonmatchesAcrossLetters.png'],'png');
% plot2svg([colorsavdir '/' group '.bootstrapciwidthforeachletter.' names{colorcat+1} '.distofnonmatchesAcrossLetters.svg'],gcf);
% close(gcf);
% 
% 
% 
% % make one averaged across all letters
%   
% figure('Name',['width of ootstrap confidence interval across letters '],...
%     'Color','White','Position',get(0,'ScreenSize'));
% 
%     hold on;
%     plot(btsamplesizes,median(squeeze(ptile(3,:,:))-squeeze(ptile(1,:,:))));
%     ylabel(['95% confidence interval on prob letter is ' names{colorcat+1}]);
%     xlabel('bootstrap sample size');
%     set(gca,'YLim',[0 .1]);
% 
% 
% saveas(gcf,[colorsavdir '/' group '.bootstrapciallletters.' names{colorcat+1} '.distofnonmatchesAcrossLetters.png'],'png');
% plot2svg([colorsavdir '/' group '.bootstrapciallletters.' names{colorcat+1} '.distofnonmatchesAcrossLetters.svg'],gcf);
% close(gcf);
% 
% 



end
















% 
% [nancor nanp] = corr(log(lewandfq),pctnansbyletter','rows','pairwise','type','spearman');
% figure;
% 
% plot(log(lewandfq),pctnansbyletter','ro');
% hold on;
% text(log(lewandfq),pctnansbyletter',letters);

