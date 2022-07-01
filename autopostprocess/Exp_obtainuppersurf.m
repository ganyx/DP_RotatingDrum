clc;clear;

addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\github_repo');%for export figure

filenum=[157:172]';
filedir=('d:\EXP RoDrtest\Exp\sharpened\');
for fi=1:size(filenum,1)
disp(['C0' num2str(filenum(fi,1)) 'ONGOIN.']);
%%
tic
load([filedir 'C0' num2str(filenum(fi,1)) '\C0' num2str(filenum(fi,1)) '_particleboundary_2.mat']);
if filenum(fi,1)<156
    load(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\' ...
    'calibration_coe_c0149-c0155.mat']);
else
    load(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\' ...
    'calibration_coe_c0157-c0172.mat']);
end
bxtmp=[];bytmp=[];
for j=1:100
bxtmp=[bxtmp;cell2mat(bx(j))];bytmp=[bytmp;cell2mat(by(j))];
end
toc
%%
tic
xtmp=[1:1:1920];
ytmp=flip([1:1:1080]');
xtmprepmat=repmat(xtmp,size(ytmp,1),1);
ytmprepmat=repmat(ytmp,1,size(xtmp,2));
xreal=xtmprepmat*pixeltorealdimension;
yreal=ytmprepmat*pixeltorealdimension;
toc
%%
tic
xi=[xreal(:) yreal(:)];
[f,ep]=ksdensity([bxtmp bytmp],xi); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),size(xreal,1),size(yreal,2));
Y = reshape(ep(:,2),size(xreal,1),size(yreal,2));
Z = reshape(f,size(xreal,1),size(yreal,2));
toc
tic
for i=1:size(Z,2)
    if sum(Z(:,i))>0
    Zindex=find(Z(:,i)==max(Z(:,i)));
    index_extreme_density(i,1)=Zindex(1);
    Xtmp=X(:,i);
    Ytmp=Y(:,i);
    surface_x(i,1)=Xtmp(index_extreme_density(i,1),1);
    surface_y(i,1)=Ytmp(index_extreme_density(i,1),1);
    else
    index_extreme_density(i,1)=nan;
    end
end
cd([filedir 'C0' num2str(filenum(fi,1))]);
figure(1)
fig=contourf(X,Y,Z,50,'LineColor','none');
hold on;
plot(surface_x,surface_y,'k.')
export_fig uppersurfaceline_contour_2.jpg -native
close figure 1
figure(2)
expfig=imread([filedir 'C0' num2str(filenum(fi,1)) '\C0'  num2str(filenum(fi,1)) '_00001.jpg']);
imshow(expfig)
hold on;plot(surface_x/pixeltorealdimension,(surface_y/pixeltorealdimension-1080)*-1,'k.');
export_fig uppersurfaceline_expfig_2.jpg -native
close figure 2
disp(['C0' num2str(filenum(fi,1)) 'DONE.']);
datafilename=['C0' num2str(filenum(fi)) '_surfacepoints_2.mat'];
save(datafilename,'surface_x');
save(datafilename,'surface_y','-append');
clear surface_x surface_y
% strr_thresh=2.5;
% [strr_bound_tmp]=contours(xreal,yreal,strrave2,[strr_thresh strr_thresh]);
% strr_bound_tmp=strr_bound_tmp';
% circle_index=find(strr_bound_tmp(:,2)==max(strr_bound_tmp(:,2)));
% strr_bound_tmp1=strr_bound_tmp(circle_index+1:circle_index+max(strr_bound_tmp(:,2)),:);
% figure(2);hold on;
% plot(strr_bound_tmp1(:,1),strr_bound_tmp1(:,2),'k.');
% % strr_bound_tmp1(circle_index+1:end,:)=[];
% % strr_bound_tmp1=strr_bound_tmp1(2:end,:);
toc
end