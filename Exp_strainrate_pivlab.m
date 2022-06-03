clc;
clear;
% plot strain rate from results of PIVlab
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function\');
filedir=(['e:\EXP RoDrtest\Exp\analysis data\']);
fileinfo=dir([filedir '*.mat']);
[fileinfo]=sortnamebysequence(fileinfo);
[Rtxt,Ctxt]=size(fileinfo);
%%
for fi=6
%     :Rtxt
    %% load data
    data = load([filedir fileinfo(fi).name]);
    filename=fileinfo(fi).name;
    expnum_tmp=split(filename,'_');
    expnum=expnum_tmp(1);
    %% extract data
    cellstrr=data.strain_rate;% strain rate cell
    cellsher=data.shear_rate;% strain rate cell
    cellx=data.x;
    celly=data.y;
    cellu=data.u_original;
    cellv=data.v_original;
    utotal=zeros(size(cell2mat(cellx(1)),1),size(cell2mat(cellx(1)),2));
    vtotal=zeros(size(cell2mat(cellx(1)),1),size(cell2mat(cellx(1)),2));
    strrtotal=zeros(size(cell2mat(cellx(1)),1),size(cell2mat(cellx(1)),2));
    shertotal=zeros(size(cell2mat(cellx(1)),1),size(cell2mat(cellx(1)),2));
    totalframes=size(cellx,1);
        for si=1:totalframes
%             if si==1
                tmpu=cell2mat(cellu(si));
                tmpu(isnan(tmpu))=0;
                tmpv=cell2mat(cellv(si));
                tmpv(isnan(tmpv))=0;
                tmpstrr=cell2mat(cellstrr(si));
                tmpsher=cell2mat(cellsher(si));
%             else
                utotal=utotal+tmpu;
                vtotal=vtotal+tmpv;
%                 strrtotal=strrtotal+tmpstrr;
%                 shertotal=shertotal+tmpsher;
%             end
        end

     uave=utotal/totalframes/1;
     vave=vtotal/totalframes/1;
%       strrave=strrtotal/totalframes*1000000;
     strrave1=strrtotal/totalframes/1;
     sherave1=shertotal/totalframes/1;
    %% ========= plot strain rate contour: accumulate after strain rate each frame ===========
%     figure(1)
% %     strrave(strrave>10000)=10000;
% %     strrave(strrave<-10000)=-10000;
%     [c1,h1]=contourf(cell2mat(cellx(si)),(cell2mat(celly(si))-max(cell2mat(celly(si)),[],'all'))*-1,sherave1,200);
%     set(h1, 'edgecolor','none');
% %     set(h,'Ydir','reverse');
    %% ========= plot strain rate contour: accumulate u&v then obtain strain rate ===========
    figure(2)
    % % ------ make coordinate same as simulation x+ left to right; y+ bottom to top-------
    xreal=flip((cell2mat(cellx(si))-max(cell2mat(cellx(si)),[],'all'))*-1,2);
%     xreal=cell2mat(cellx(si));
    yreal=flip((cell2mat(celly(si))-max(cell2mat(celly(si)),[],'all'))*-1,2);
%     yreal=cell2mat(celly(si));
%     Vxavereal=-flip(uave,2);
%     Vyavereal=flip(vave,2);
    Vxavereal=uave;
    Vyavereal=-vave;
    
    strrave2=strain(xreal,yreal,Vxavereal,Vyavereal);
    [c2,h2]=contourf(xreal,yreal,strrave2,200);
    set(h2, 'edgecolor','none');
    axis equal;
    xlim([0.05 0.3]);
    colorbar;
end
%%
load(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\' char(expnum) '_particleboundary.mat']);
bxtmp=[];bytmp=[];
for j=1:100
bxtmp=[bxtmp;cell2mat(bx(j))];bytmp=[bytmp;cell2mat(by(j))];
end
xi=[xreal(:) yreal(:)];
[f,ep]=ksdensity([bxtmp bytmp],xi); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),size(xreal,1),size(yreal,2));
Y = reshape(ep(:,2),size(xreal,1),size(yreal,2));
Z = reshape(f,size(xreal,1),size(yreal,2));
figure(3)
fig=contourf(X,Y,Z,50,'LineColor','none');
hold on;

for i=1:size(Z,2)
    if sum(Z(:,i))>0
    index_extreme_density(i,1)=find(Z(:,i)==max(Z(:,i)));
    Xtmp=X(:,i);
    Ytmp=Y(:,i);
    surface_x(i,1)=Xtmp(index_extreme_density(i,1),1);
    surface_y(i,1)=Ytmp(index_extreme_density(i,1),1);
    else
    index_extreme_density(i,1)=nan;
    end
end
plot(surface_x,surface_y,'k.');

% figure(4)
[strr_bound_tmp]=contours(xreal,yreal,strrave2,[5 5]);
% plot(x_extremedensity,y_extremedensity,'k.');
strr_bound_tmp=strr_bound_tmp';
circle_index=find(strr_bound_tmp(:,2)==max(strr_bound_tmp(:,2)));
strr_bound_tmp1=strr_bound_tmp(circle_index+1:circle_index+max(strr_bound_tmp(:,2)),:);
% strr_bound_tmp1(circle_index+1:end,:)=[];
% strr_bound_tmp1=strr_bound_tmp1(2:end,:);
%% method 1 extract shear region lower boundary: particle to surface distance
% cali=0.236/1295.64;
% pixeltorealdimension=cali/0.96;
% maxdist_index=[];
% dist_to_surfacecell=[];
% count=0;
% for sbi=1:size(surface_x,1)
%     strr_index=find(strr_bound_tmp1(2:end,1)==surface_x(sbi));
%     if ~isempty(strr_index)
%         count=count+1;
%         strr_xtmp=strr_bound_tmp1(strr_index,1);
%         strr_ytmp=strr_bound_tmp1(strr_index,2);
%         dist_to_surface=strr_ytmp-surface_y(sbi);
%         if size(dist_to_surface,1)<2
%             count=count-1;
%         else
%             dist_to_surfacecell{count,1}=strr_ytmp-surface_y(sbi);
%             maxdist_index(count,1)=strr_index(find(dist_to_surface==min(dist_to_surface)));
%         end
%     else
%         continue
%     end
% 
% end
% strr_bound=strr_bound_tmp(maxdist_index,:);
% boundary_in_pixel=[surface_x,surface_y;strr_bound]/pixeltorealdimension;

%% method 2 extract shear region lower boundary: clusterdata function
cali=0.236/1295.64;
pixeltorealdimension=cali/0.96;
dp=0.002;%particle size
xws=2*dp;%x window size
steps=1*dp;%moving step size
[strr_botbound]=strainrate_bottombound(strr_bound_tmp1,xws,steps);
boundary_in_pixel=[surface_x,surface_y;strr_botbound]/pixeltorealdimension;
boundary_in_pixel(:,2)=(boundary_in_pixel(:,2)-1080)*-1;
datafilename=[char(expnum) 'strinrate_boundary.mat'];
save(datafilename,'boundary_in_pixel');

%% fitting and find the strainrate depth
marker_color_type='ko';
dpfit=0.002;%particle size
xwsfit=10*dpfit;%x window size
stepsfit=1*dpfit;%moving step size
[strr_depth,strr_depth_normbydp,maxdepth_norm,swfit]=strr_depth_fun(strr_botbound,surface_x,surface_y,dp,xwsfit,stepsfit);
datafilename=([char(expnum) '_strainrate_depth.mat']);
save(datafilename,'strr_depth_normbydp');
save(datafilename,'maxdepth_norm','-append');
save(datafilename,'swfit','-append');

%%
function out=strain(x,y,vx,vy)
hx = x(1,:);
hy = y(:,1);
[pvxx, pvxy] = gradient(vx, hx, hy);
[pvyx, pvyy] = gradient(vy, hx, hy); %#ok<*ASGLU>
out = 0.5*(-pvxy+pvyx);
end

function out=shear(x,y,u,v)
hx = x(1,:);
hy = y(:,1);
[junk, py] = gradient(u, hx, hy);
[qx, junk] = gradient(v, hx, hy);
out= 0.5*(qx+py);
end



function [strr_bound]=strainrate_bottombound(strr_bound_tmp1,xws,steps)
xrange=[min(strr_bound_tmp1(:,1)) max(strr_bound_tmp1(:,1))];%fron extreme left point to right point
sw=[xrange(1):steps:xrange(2)]';%sample window, set steps to xws to extract boundary points (otherwise overlap points)
strr_bound=[];
for wi=1:size(sw,1)
    window=[sw(wi)-xws/2 sw(wi)+xws/2];
    points_index1=find(strr_bound_tmp1(:,1)>window(1));
    points_index2=find(strr_bound_tmp1(:,1)<window(2));
    points_index=points_index1(ismember(points_index1,points_index2));
    points_tmp=strr_bound_tmp1(points_index,:);
    points_tmp_cell{wi,1}=points_tmp;
    cluster_index=clusterdata(points_tmp,2);
    ycluster1=mean(points_tmp(cluster_index==1,2));
    ycluster2=mean(points_tmp(cluster_index==2,2));
    current_strr_bot=[];
    if ycluster1>ycluster2
        strr_bound=[strr_bound;points_tmp(cluster_index==2,:)];
    else
        strr_bound=[strr_bound;points_tmp(cluster_index==1,:)];
    end
end

end

function [strr_depth,strr_depth_normbydp,maxdepth_norm,swfit]=strr_depth_fun(strr_botbound,surface_x,surface_y,dp,xwsfit,stepsfit)
xrangefit=[min(strr_botbound(:,1)) max(strr_botbound(:,1))];%fron extreme left point to right point
swfit=[xrangefit(1):stepsfit:xrangefit(2)]';%sample window, set steps to xws to extract boundary points (otherwise overlap points)
fitangle_a=[];
fitangle_b=[];
gof_fitangle=[];
fitdatamidx=[];
fitdatamidy=[];
% %--------------angle fitting settings----------------
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.395894484205228 0.936450101590675];
% %--------------end angle fitting settings----------------

for fiti=1:size(swfit,1)
%     :size(swfit,1)
    % Set up fittype and options.
    windowfit=[swfit(fiti)-xwsfit/2 swfit(fiti)+xwsfit/2];
    points_index1=find(strr_botbound(:,1)>windowfit(1));
    points_index2=find(strr_botbound(:,1)<windowfit(2));
    points_index=points_index1(ismember(points_index1,points_index2));
    points_tmp=strr_botbound(points_index,:);
        [fitresult, gof] = fit(points_tmp(:,1),points_tmp(:,2),ft,opts);
        fitangle_a(fiti,1)=atan(fitresult.a);%angle in radian
        fitangle_b(fiti,1)=fitresult.b;%intersection between y axis
        gof_fitangle(fiti,1)=gof.rsquare;
%         DAOR(samplei,1)=max(fitangle_a);     
%         DAOR_gof(samplei,1)=gof_fitangle(fitangle_a==max(fitangle_a));
    %% mid point of bottom
    midbot=[swfit(fiti) mean(points_tmp(:,2))];
    %% left and right point based on bottom mid point
    botleft_coor=[midbot(1)-xwsfit/2 midbot(2)+xwsfit/2*fitangle_a(fiti,1)];
    botright_coor=[midbot(1)+xwsfit/2 midbot(2)+xwsfit/2*fitangle_a(fiti,1)];
    %% left and right line and second point on each line
    lr_tan=tan(fitangle_a(fiti,1)+pi/2);
    left_b=botleft_coor(2)-botleft_coor(1)*lr_tan;
    right_b=botright_coor(2)-botright_coor(1)*lr_tan;
    left_2nd_coor=[0 left_b];
    right_2nd_coor=[0 right_b];
    %% dist to lines of surface points & extract between lines
    dist_linetoline=sqrt(sum((botleft_coor-botright_coor).^2));
    dist_surf_left_right(:,1)=point_to_line_distance([surface_x surface_y],botleft_coor,left_2nd_coor);
    dist_surf_left_right(:,2)=point_to_line_distance([surface_x surface_y],botright_coor,right_2nd_coor);
    betweenlines_indexleft=find(dist_surf_left_right(:,1)<dist_linetoline);
    betweenlines_indexright=find(dist_surf_left_right(:,2)<dist_linetoline);
    betweenlines_index=betweenlines_indexleft(ismember(betweenlines_indexleft,betweenlines_indexright));
    betweenlines_pointscoor=[surface_x(betweenlines_index) surface_y(betweenlines_index)];
    midtop=[mean(betweenlines_pointscoor(:,1)) mean(betweenlines_pointscoor(:,2))];
    %% flow region depth; dist between midbot and midtop
    strr_depth(fiti,1)=sqrt(sum((midtop-midbot).^2));
    strr_depth_normbydp=strr_depth/dp;
end
maxdepth_norm=max(strr_depth_normbydp);
% figure(1);hold on;plot(swfit,fitangle_a/pi*180,marker_color_type);
% figure(3);hold on;plot(swfit,strr_depth_normbydp,marker_color_type);
end
