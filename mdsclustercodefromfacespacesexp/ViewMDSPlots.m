% makes plots organized by view.  essentially takes the output of
% SsRateAllPairsCompareSubsets and generates the same plots except that the
% former assumes there are 6 groups (ids) which cluster by view, and this
% output assumes there are 4 groups (view) which cluster the ids.  will
% only work at the moment if SsRateAllPairsCompareSubsets is run first

% need to resort faces{i} and ids{i} so that they are organized by view
% instead of identity

% since the number of view we want to analyze is variable will have to
% derive it from variable compare which says which views to analyze
whitebg([1 1 1]);



v_ids=[];
vn_labels = cell(length(new_labels),1);
v_faces = faces;
ind=1;
for i=1:length(compare) %for each view
    for j=1:6 %grab from each identity
        v_ids =[v_ids; ids((i+(j-1)*length(compare)),:)]; %reindex image names
        vn_labels{ind} = new_labels{v_ids(ind,2)};%reindex image labels, these are the ones supplied by altnames
        v_faces{ind} = faces{i+(j-1)*length(compare)}; %reindex images
        ind=ind+1;
    end
end



% want to plot arbitrary pair of dimensions in a figure
% need to know dimensionality of the solution
% and which two dimensions to plot

num_dims_insol = 2;
dim1 = 1;
dim2 = 2;


%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

% fill y dim of plot if solution is only 1 dimensional
if num_dims_insol == 1
    sol_to_plot(:,2)=0;
end


% resort solutions to match new list order

v_sol_to_plot = [];

% for i=1:4
for i=1:length(compare)
    for j=1:6
        v_sol_to_plot =[v_sol_to_plot; sol_to_plot((i+(j-1)*length(compare)),:)]; %reindex solutions
    end
end


%now plot chosen dimensions



imsize = 0.2;
fig_posn = [0 0 1000 1000];

n=n+1;
F=figure(n); 


%plot first face
image([v_sol_to_plot(1,dim1)-imsize;v_sol_to_plot(1,dim1)+imsize],[v_sol_to_plot(1,dim2)-imsize; v_sol_to_plot(1,dim2)+imsize],v_faces{1});
hold on;

%now the rest  need to paramterize this so that it only uses as many faces
%as in subset
for ctr = 2:length(ids)
   image([v_sol_to_plot(ctr,dim1)-imsize; v_sol_to_plot(ctr,dim1)+imsize],[v_sol_to_plot(ctr,dim2)-imsize; v_sol_to_plot(ctr,dim2)+imsize],v_faces{ctr}); 
end


% number of views
nv = length(compare);
ni = 6;  %number of ids

% colors for plot, will only use first 4 or less
labelmap = [1 0 0; 0 .4 0; .5 0 .75; 0 0 1; 0 1 0; 0 .5 0];

% want to loop nv times giving ni colors in a row
% for i = 1:nv
% %     solutions should go 1-6,7-12 nv times
%     text(v_sol_to_plot(i*ni-(ni-1):i*ni,dim1),v_sol_to_plot(i*ni-(ni-1):i*ni,dim2),num2str(v_ids(i*ni-(ni-1):i*ni)'),'Color',labelmap(i,:), 'FontSize',18);
% end

for i = 1:nv
% %     solutions should go 1-6,7-12 nv times
    text(v_sol_to_plot(i*ni-(ni-1):i*ni,dim1),v_sol_to_plot(i*ni-(ni-1):i*ni,dim2),vn_labels(i*ni-(ni-1):i*ni)','Color',labelmap(i,:), 'FontSize',30,'FontWeight','Bold');
end




xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,0,1],'FontSize',18);
title(['mds on pixel differences: '    num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);

hold off;
%find how big the plot needs to be
min_x = min(v_sol_to_plot(:,dim1));
max_x = max(v_sol_to_plot(:,dim1));
min_y = min(v_sol_to_plot(:,dim2));
max_y = max(v_sol_to_plot(:,dim2));

vals = abs([min_x max_x min_y max_y]);

%largest distance + some wiggle room
val = max(vals) + imsize;
%set axes
axis([-val val -val val]);
colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% WANT TO PLOT ARBITRARY PAIR OF DIMENSIONS IN A 3DFIGURE
% NEED TO KNOW DIMENSIONALITY OF THE SOLUTION
% AND WHICH TWO DIMENSIONS TO PLOT

num_dims_insol = 3;
dim1 = 1;
dim2 = 2;
dim3 = 3;
%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

% resort solutions to match new list order

v_sol_to_plot = [];

for i=1:4
    for j=1:6
        v_sol_to_plot =[v_sol_to_plot; sol_to_plot((i+(j-1)*4),:)]; %reindex solutions
    end
end




%now plot chosen dimensions



imsize = 0.2;
fig_posn = [0 0 1000 1000];

n=n+1;
F=figure(n); 


%plot first face
image([v_sol_to_plot(1,dim1)-imsize;v_sol_to_plot(1,dim1)+imsize;],[v_sol_to_plot(1,dim2)-imsize; v_sol_to_plot(1,dim2)+imsize],v_faces{1});
hold on;

%now the rest
for ctr = 2:length(ids)
   image([v_sol_to_plot(ctr,dim1)-imsize; v_sol_to_plot(ctr,dim1)+imsize],[v_sol_to_plot(ctr,dim2)-imsize; v_sol_to_plot(ctr,dim2)+imsize],v_faces{ctr}); 
end




for i = 1:nv
    text(v_sol_to_plot(i*ni-(ni-1):i*ni,dim1),v_sol_to_plot(i*ni-(ni-1):i*ni,dim2),v_sol_to_plot(i*ni-(ni-1):i*ni,dim3) ,num2str(v_ids(i*ni-(ni-1):i*ni)'),'Color',labelmap(i,:), 'FontSize',18);
end



xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,1,0],'FontSize',18);
 zlabel(['dimension' num2str(dim3)],'Color',[0,1,0],'FontSize',18);
 
title(['mds on ' subname ' ratings : '   num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);


%find how big the plot needs to be
min_x = min(v_sol_to_plot(:,dim1));
max_x = max(v_sol_to_plot(:,dim1));
min_y = min(v_sol_to_plot(:,dim2));
max_y = max(v_sol_to_plot(:,dim2));
min_z = min(v_sol_to_plot(:,dim2));
max_z = min(v_sol_to_plot(:,dim3));

vals = abs([min_x max_x min_y max_y min_z max_z]);

%largest distance + some wiggle room
val = max(vals) + imsize+1;
%set axes
axis([-val val -val val -val val]);
colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);

hold off;




cd (thepath);



%plot points surrounded by convex hull.  should allow visualization of
%groupings more easily

% first need num dimensions in solution and dimensions to plot
n=n+1;
figure(n);


num_dims_insol = 3;
dim1 = 1;
dim2 = 2;
dim3 = 3;
%get solution of right dimensionality
sol_to_plot = dim_solutions{num_dims_insol};

% resort solutions to match new list order

v_sol_to_plot = [];

for i=1:4
    for j=1:6
        v_sol_to_plot =[v_sol_to_plot; sol_to_plot((i+(j-1)*4),:)]; %reindex solutions
    end
end




% set up colormap for figure to get colors right

colormap([1 0 0
        .5 0 0
        0 0 1
        0 0 .5
        0 1 0
        0 .5 0]);

    caxis([1,6]);

col = [1 1 1 1
        2 2 2 2
        3 3 3 3 
        4 4 4 4
        5 5 5 5
        6 6 6 6];

%now find convex hull for each of the identities



for i = 1:nv
     patch(v_sol_to_plot(i*ni-(ni-1):i*ni,1),v_sol_to_plot(i*ni-(ni-1):i*ni,2),v_sol_to_plot(i*ni-(ni-1):i*ni,3),labelmap(i,:));alpha(.5);
%     trisurf(hull,v_sol_to_plot(i*ni-(ni-1):i*ni,1),v_sol_to_plot(i*ni-(ni-1):i*ni,2),v_sol_to_plot(i*ni-(ni-1):i*ni,3),col(i,:)); alpha(.5);
end


hold on;


xlabel(['dimension' num2str(dim1)],'Color',[1,0,0],'FontSize',18);  ylabel(['dimension' num2str(dim2)],'Color',[0,1,0],'FontSize',18);
 zlabel(['dimension' num2str(dim3)],'Color',[0,1,0],'FontSize',18);
 
title(['mds on ' subname ' ratings : '   num2str(num_dims_insol) ' dimensional solution'],'FontSize',24);

%find how big the plot needs to be
min_x = min(v_sol_to_plot(:,dim1));
max_x = max(v_sol_to_plot(:,dim1));
min_y = min(v_sol_to_plot(:,dim2));
max_y = max(v_sol_to_plot(:,dim2));
min_z = min(v_sol_to_plot(:,dim2));
max_z = min(v_sol_to_plot(:,dim3));

vals = abs([min_x max_x min_y max_y min_z max_z]);

%largest distance + some wiggle room
val = max(vals) + imsize+1;
%set axes
axis([-val val -val val -val val]);
% colormap([[0:1/255:1]', [0:1/255:1]',[0:1/255:1]'])
set(F,'Position',fig_posn);

hold off;


colormap('default');




% make a shepard diagram  plotting proximities in the solution against distances in the data 









