% Multidimensional scaling for color matching data.  At the moment runs on
% the average CIEL*a*b* distance across subjects



% Get scaling solution for mean distance matrix.  
% number of iterations to do before giving up
opts = statset('MaxIter',2000)
% max number of dimensions to solve for
maxdms = 25;

% clear stress, mean_mds_solution, disparities;


criteria = {'stress', 'sstress', 'metricstress','metricsstress','sammon','strain'};

% for each type of criterion
for c=1:length(criteria)
    
%     solve for all dimensionalities and plot shepard plot
    figure('Name',['Shepard Plots from MDS: ' criteria{c}],'Color',[1 1 1],'Position',get(0,'ScreenSize'));
    
    for dms=1:maxdms
        [mean_mds_solution{dms},stress(dms),disparities{dms}] ...
            = mdscale(squareform(meanlabdist),...
            dms,'Criterion',criteria{c},'Start','random','Options',opts);
        
        % make a shepard plot
        subplot(5,5,dms);
        
        plot(meanlabdist,pdist(mean_mds_solution{dms}),'ro','MarkerFaceColor',[1 0 0],'MarkerSize',5');
        xlabel('original distances');
        ylabel('solution distances');
        hold on;
        box off;
        %     add x=y line
        plot(0:1:max(meanlabdist),0:1:max(meanlabdist),'k--');
        axis equal;
        title([num2str(dms) ' dimensions']);
        
    end
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %
    % plot stress x dimensionality
    figure('Name',['Scree Plots: ' criteria{c}],'Color',[1 1 1]);
    
    plot(1:length(stress),stress,'bo','MarkerSize',10,'MarkerFaceColor',[0 0 1]);
    xlabel('number of dimensions in solution','FontSize',18);
    ylabel('stress','FontSize',18);
    title('goodness of fit for solutions of different dimensionality');
    box off;
end



% only the sammon metric gets under .1 stress in 3 dimensions. the others
% require more.  in all cases the shepard plots look funky until you have
% 15 or more dimensions really, tending to underestimate short distances
% and overestimate long ones.  


% let's just get the solutions with the sammon distance metric

for dms=1:maxdms
    [mean_mds_solution{dms},stress(dms),disparities{dms}] ...
        = mdscale(squareform(meanlabdist),...
        dms,'Criterion','sammon','Start','random','Options',opts);
    
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 


dims2plot = [1 2 3];
diminsol = 15;

figure('Name','MDS Solution','Color',[1 1 1],'Position',get(0,'ScreenSize'));

scatter3(mean_mds_solution{diminsol}(:,dims2plot(1)),...
    mean_mds_solution{diminsol}(:,dims2plot(2)),...
    mean_mds_solution{diminsol}(:,dims2plot(3)),'r');

hold on;

text(mean_mds_solution{diminsol}(:,dims2plot(1)),...
    mean_mds_solution{diminsol}(:,dims2plot(2)),...
    mean_mds_solution{diminsol}(:,dims2plot(3)), letters,'FontSize',14);

xlabel('dim1');
ylabel('dim2');
zlabel('dim3');




% let's do it as 3 subplots
dims2plot = [1 2 3];
dimsinsol = 3;

figure('Name',['MDS Solution paired dimensions ' num2str(dimsinsol) ' dimensional solution'],...
    'Color',[1 1 1],'Position',get(0,'ScreenSize'));

subplot(1,3,1);

plot(mean_mds_solution{diminsol}(:,dims2plot(1)),...
    mean_mds_solution{diminsol}(:,dims2plot(2)),'ro','MarkerFaceColor','r',...
    'MarkerSize',10);

hold on;

text(mean_mds_solution{diminsol}(:,dims2plot(1)),...
    mean_mds_solution{diminsol}(:,dims2plot(2)),...
    letters,'FontSize',14);

box off
% axis square
% 
axis equal

% axis image

xlabel(['dimension ' num2str(dims2plot(1))]);
ylabel(['dimension ' num2str(dims2plot(2))]);

% % % % %

subplot(1,3,2);

plot(mean_mds_solution{diminsol}(:,dims2plot(1)),...
    mean_mds_solution{diminsol}(:,dims2plot(3)),'ro','MarkerFaceColor','r',...
    'MarkerSize',10);

hold on;

text(mean_mds_solution{diminsol}(:,dims2plot(1)),...
    mean_mds_solution{diminsol}(:,dims2plot(3)),...
    letters,'FontSize',14);
box off
axis equal
% axis square

xlabel(['dimension ' num2str(dims2plot(1))]);
ylabel(['dimension ' num2str(dims2plot(3))]);


% % % % % %
subplot(1,3,3);

plot(mean_mds_solution{diminsol}(:,dims2plot(2)),...
    mean_mds_solution{diminsol}(:,dims2plot(3)),'ro','MarkerFaceColor','r',...
    'MarkerSize',10);
hold on;

text(mean_mds_solution{diminsol}(:,dims2plot(2)),...
    mean_mds_solution{diminsol}(:,dims2plot(3)),...
    letters,'FontSize',14);
box off
axis equal
% axis square

xlabel(['dimension ' num2str(dims2plot(2))]);
ylabel(['dimension ' num2str(dims2plot(3))]);












% %get solution of right dimensionality
% sol_to_plot = dim_solutions{num_dims_insol};
% 
% % fill y dim of plot if solution is only 1 dimensional
% if num_dims_insol == 1
%     sol_to_plot(:,2)=0;
% end
% 
% 
% %now plot chosen dimensions
% 
% 
% 
% imsize = 0.3;
% fig_posn = [0 0 1000 1000];
% 
% n=n+1;
% F=figure(n); 
% 
% 
% %plot first face
% image([sol_to_plot(1,dim1)-imsize;sol_to_plot(1,dim1)+imsize],[sol_to_plot(1,dim2)-imsize; sol_to_plot(1,dim2)+imsize],faces{1});
% hold on;
% 
% %now the rest  need to paramterize this so that it only uses as many faces
% %as in subset
% for ctr = 2:length(ids)
%    image([sol_to_plot(ctr,dim1)-imsize; sol_to_plot(ctr,dim1)+imsize],[sol_to_plot(ctr,dim2)-imsize; sol_to_plot(ctr,dim2)+imsize],faces{ctr}); 
% end
% 
% 
% % will always be 6 identities
% nv = length(compare);
% 
% % colors for plot
% labelmap = [1 0 0; .5 0 .5; 0 0 1; 1 .5 0; 1 1 0; 0 1 0];
% 
% 
% % for i = 1:length(ids)/nv
% %     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2),num2str([i*nv-(nv-1):i*nv]'),'Color',labelmap(i,:), 'FontSize',18);
% % end
% 
% 
% 
% % for i = 1:length(ids)/nv
% %     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2) ,num2str(ids(i*nv-(nv-1):i*nv)'),'Color',labelmap(i,:), 'FontSize',24,'FontWeight','Bold');
% % end
% 
% 
% % alternatively plot names from morph space (supplied by variable altname
% % that would be loaded with the data
% for i = 1:length(ids)/nv
%     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2) ,new_labels(i*nv-(nv-1):i*nv)','Color',labelmap(i,:), 'FontSize',20,'FontWeight','Bold');
% end
% 
% xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,1,0],'FontSize',18);
% title(['mds on ' subname ' ratings : '    num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);
% 
% %find how big the plot needs to be
% min_x = min(sol_to_plot(:,dim1));
% max_x = max(sol_to_plot(:,dim1));
% min_y = min(sol_to_plot(:,dim2));
% max_y = max(sol_to_plot(:,dim2));
% 
% vals = abs([min_x max_x min_y max_y]);
% 
% %largest distance + some wiggle room
% val = max(vals) + imsize;
% %set axes
% axis([-val val -val val]);
% colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
% set(F,'Position',fig_posn);
% hold off;
% 
% 
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % WANT TO PLOT ARBITRARY PAIR OF DIMENSIONS IN A 3DFIGURE
% % NEED TO KNOW DIMENSIONALITY OF THE SOLUTION
% % AND WHICH TWO DIMENSIONS TO PLOT
% 
% num_dims_insol = 3;
% dim1 = 1;
% dim2 = 2;
% dim3 = 3;
% %get solution of right dimensionality
% sol_to_plot = dim_solutions{num_dims_insol};
% 
% %now plot chosen dimensions
% 
% 
% 
% imsize = 0.3;
% fig_posn = [0 0 1000 1000];
% 
% n=n+1;
% F=figure(n); 
% 
% 
% %plot first face
% image([sol_to_plot(1,dim1)-imsize;sol_to_plot(1,dim1)+imsize;],[sol_to_plot(1,dim2)-imsize; sol_to_plot(1,dim2)+imsize],faces{1});
% hold on;
% 
% %now the rest
% for ctr = 2:length(ids)
%    image([sol_to_plot(ctr,dim1)-imsize; sol_to_plot(ctr,dim1)+imsize],[sol_to_plot(ctr,dim2)-imsize; sol_to_plot(ctr,dim2)+imsize],faces{ctr}); 
% end
% 
% 
% 
% % plot names of images on screen
% % for i = 1:length(ids)/nv
% %     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2),sol_to_plot(i*nv-(nv-1):i*nv,dim3) ,num2str(ids(i*nv-(nv-1):i*nv)'),'Color',labelmap(i,:), 'FontSize',18);
% % end
% 
% % alternatively plot names from morph space (supplied by variable altname
% % that would be loaded with the data
% for i = 1:length(ids)/nv
%     text(sol_to_plot(i*nv-(nv-1):i*nv,dim1),sol_to_plot(i*nv-(nv-1):i*nv,dim2),sol_to_plot(i*nv-(nv-1):i*nv,dim3) ,new_labels(i*nv-(nv-1):i*nv)','Color',labelmap(i,:), 'FontSize',20,'FontWeight','Bold');
% end
% 
% % label axes
% xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,0,1],'FontSize',18);
%  zlabel(['dimension' num2str(dim3)],'Color',[0,1,0],'FontSize',18);
%  
% title(['mds on ' subname ' ratings : '   num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);
% 
% %find how big the plot needs to be
% min_x = min(sol_to_plot(:,dim1));
% max_x = max(sol_to_plot(:,dim1));
% min_y = min(sol_to_plot(:,dim2));
% max_y = max(sol_to_plot(:,dim2));
% min_z = min(sol_to_plot(:,dim2));
% max_z = min(sol_to_plot(:,dim3));
% 
% vals = abs([min_x max_x min_y max_y min_z max_z]);
% 
% %largest distance + some wiggle room
% val = max(vals) + imsize+1;
% %set axes
% axis([-val val -val val -val val]);
% colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
% set(F,'Position',fig_posn);
% hold off;
% 
% 
% 
% 
% 
% cd (thepath);
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% 
% %plot points surrounded by convex hull.  should allow visualization of
% %groupings more easily
% 
% % first need num dimensions in solution and dimensions to plot
% n=n+1;
% figure(n);
% 
% 
% num_dims_insol = 3;
% dim1 = 1;
% dim2 = 2;
% dim3 = 3;
% %get solution of right dimensionality
% sol_to_plot = dim_solutions{num_dims_insol};
% 
% % set up colormap for figure to get colors right
% 
% colormap([1 0 0; .5 0 .5; 0 0 1; 1 .5 0; 1 1 0; 0 1 0]);
% 
%     caxis([1,6]);   
% 
% col = [1 1 1 1
%         2 2 2 2
%         3 3 3 3 
%         4 4 4 4
%         5 5 5 5
%         6 6 6 6];
% 
% %now find convex hull for each of the identities
% 
% % 
% % hull=[sol_to_plot(:,1),sol_to_plot(:,2),sol_to_plot(:,3)];
% % for i = 1:length(ids)/nv
% %      patch(sol_to_plot(i*nv-(nv-1):i*nv,1),sol_to_plot(i*nv-(nv-1):i*nv,2),sol_to_plot(i*nv-(nv-1):i*nv,3),labelmap(i,:));alpha(.5);
% % %     trisurf(hull,sol_to_plot(i*nv-(nv-1):i*nv,1),sol_to_plot(i*nv-(nv-1):i*nv,2),sol_to_plot(i*nv-(nv-1):i*nv,3),col(i,:)); alpha(.5);
% % end
% 
% hull1 = convhulln(sol_to_plot(1:4,:));
% trisurf(hull1,sol_to_plot(1:4,1),sol_to_plot(1:4,2),sol_to_plot(1:4,3),col(1,:)); alpha(.5);
% % fill3(hull1(:,1),hull1(:,2),hull1(:,3),col1');
% hold on;
% % patch(hull1(2,1),hull1(2,2),hull1(2,3),[1 0 0]);
% % patch(hull1(3,1),hull1(3,2),hull1(3,3),[1 0 0]);
% % patch(hull1(4,1),hull1(4,2),hull1(4,3),[1 0 0]);
% 
% 
% 
% hull2 = convhulln(sol_to_plot(5:8,:));
% trisurf(hull2,sol_to_plot(5:8,1),sol_to_plot(5:8,2),sol_to_plot(5:8,3),col(2,:));alpha(.5);
% hull3 = convhulln(sol_to_plot(9:12,:));
% trisurf(hull3,sol_to_plot(9:12,1),sol_to_plot(9:12,2),sol_to_plot(9:12,3),col(3,:)); alpha(.5);
% hull4 = convhulln(sol_to_plot(13:16,:));
% trisurf(hull4,sol_to_plot(13:16,1),sol_to_plot(13:16,2),sol_to_plot(13:16,3),col(4,:)); alpha(.5);
% hull5 = convhulln(sol_to_plot(17:20,:));
% trisurf(hull5,sol_to_plot(17:20,1),sol_to_plot(17:20,2),sol_to_plot(17:20,3),col(5,:));alpha(.5);
% hull6 = convhulln(sol_to_plot(21:24,:));
% trisurf(hull6,sol_to_plot(21:24,1),sol_to_plot(21:24,2),sol_to_plot(21:24,3),col(6,:)); alpha(.5);
% hold on;
% 
% 
% xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,1,0],'FontSize',18);
%  zlabel(['dimension' num2str(dim3)],'Color',[0,1,0],'FontSize',18);
%  
% title(['mds on ' subname ' ratings : '   num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);
% 
% %find how big the plot needs to be
% min_x = min(sol_to_plot(:,dim1));
% max_x = max(sol_to_plot(:,dim1));
% min_y = min(sol_to_plot(:,dim2));
% max_y = max(sol_to_plot(:,dim2));
% min_z = min(sol_to_plot(:,dim2));
% max_z = min(sol_to_plot(:,dim3));
% 
% vals = abs([min_x max_x min_y max_y min_z max_z]);
% 
% %largest distance + some wiggle room
% val = max(vals) + imsize+1;
% %set axes
% axis([-val val -val val -val val]);
% % colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
% set(F,'Position',fig_posn);
% 
% hold off;
% 
% 
% colormap('default');
% 
% 
% 
% 
% 
% 
% 
