
% index to subjects with NaNs
subswnans = [];
% how many NaNs for each subject
numnans = [];
% for each subject
for i=1:6588
%     check for a nan
    y=sum(isnan(u_rgb(i,2,:)));
    if y>0
%     fprintf('subject %g has %g nans \n',i,y);
    subswnans = [subswnans i];
    numnans = [numnans y];
    else
        numnans = [numnans 0];
    end
end
    
% % so less than 5% of the matches are NaNs
% 6588*26=171288
%       
% sum(numnans)
% ans =
%         7352
% 7352/171288
% ans =
%      0.042922
% which is not too bad, but we need to account for them.

%  these are distributed throughout 35% of the population (which is more
%  than I would have expected)
% size(subswnans,2)/6588 =0.352

% if we do nothing, the nans get assigned the color black (0,0,0);

% the matches also get the label 0 which corresponds also to the color
% black

% however if we don't allow the interpolation at the label assignment
% stage, then things should fix themselves in the histograms.  this fix is
% in the function fpRGB2ColorsJW.m.   All the functions appear to work now.

% let's make some plots so we can see where the problems are.  
% number of histogram of subjects with n NaNs

figure('Name','histogram of subjects with n letters with no color','Color',[1 1 1]);
hist(numnans/sum(numnans),[0:26]);
xlabel('number of letters with no match');
ylabel('number of subjects');
box off;

% do this as percentages
h=hist(numnans,[0:26]);

pctnans = h/sum(h);

% 
%             0        0.648
%             1      0.12614
%             2      0.05935
%             3      0.03901
%             4     0.034153
%             5     0.025956
%             6     0.022769
%             7     0.018367
%             8     0.016697
%             9    0.0069824
%            10    0.0022769
%            11   0.00030358
%            12            0
%            13            0
%            14            0
%            15            0
%            16            0
%            17            0
%            18            0
%            19            0
%            20            0
%            21            0
%            22            0
%            23            0
%            24            0
%            25            0
%            26            0



figure('Name','percent of subjects with n letters with no color','Color',[1 1 1]);
bar(pctnans);

xlabel('number of letters with no match');
ylabel('percent of subjects');
box off;





% how is this behavior distributed across letters?
figure('Name','distribution of nans across letters','Color',[1 1 1]);

subplot(1,2,1);
% figure where nans are black and color matches are white
% get one rgb column of data for subs x letters
nansmatrix = squeeze(u_rgb(:,2,:));

nansmatrix(find(~isnan(nansmatrix)))=0;
nansmatrix(find(isnan(nansmatrix)))=1;
imagesc(nansmatrix);  
colormap(bone);
box off;
xlabel('letters');
ylabel('subjects');
set(gca,'XTick',[1:26],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);


% let's look at the distribution across letters
subplot(1,2,2);

bar(sum(nansmatrix)/6588);;
box off;
xlabel('letters');
ylabel('number of times not matched');
set(gca,'XTick',[1:26],'XTickLabel',letters);
% set(gca,'XTickLabel',letters);

