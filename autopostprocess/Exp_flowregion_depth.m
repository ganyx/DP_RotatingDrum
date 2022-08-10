clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('e:\EXP RoDrtest\Exp\sharpened\');
[exprec_num,exprec_txt,exprec_raw]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',9);
omega_xls=exprec_num(1:6,1);
watercontent=exprec_num(1,1);
dp=0.002;%particle diameter
HFR=0.002;%500 frames/second
filenum=[157,158,161]';
% filenum=[149:155,157:174]';
%%
for fi=1:size(filenum,1)
    disp(['C0' num2str(filenum(fi))]);
% figureinfo=dir([read_filedir '\*.jpg']);
%% extract calibration and Container center info
maskdir='C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\';
if filenum(fi,1)<156
    CCR=load([maskdir 'C0149-C0155_CR.mat']);
    load([maskdir 'calibration_coe_c0149-c0155.mat']);
elseif filenum(fi,1)>155 && filenum(fi,1)<177
    CCR=load([maskdir 'C0157-C0172_CR.mat']);
    load([maskdir 'calibration_coe_c0157-c0172.mat']);
elseif filenum(fi,1)>176
    CCR=load([maskdir 'C0177-C0182_CR.mat']);
    load([maskdir 'calibration_coe_c0177-c0182.mat']);
end
CCreal=CCR.CC*pixeltorealdimension;
CCreal(1,2)=(CCreal(1,2)-1080*pixeltorealdimension)*-1;

%% load upper surface data surface_x surface_y
load([filedir 'C0' num2str(filenum(fi,1)) '\C0' num2str(filenum(fi,1)) '_surfacepoints_2.mat']);
surface_x_real=surface_x-CCreal(1);
surface_y_real=surface_y-CCreal(2);

%% extrac surface points inside container
surf_disttoCC=sqrt((surface_x_real).^2+(surface_y_real).^2);
surfinsideCircle_index=find(surf_disttoCC<0.07);%0.1 is the real container radius
sicx=surface_x_real(surfinsideCircle_index);%surface in circle
sicy=surface_y_real(surfinsideCircle_index);

%% Locate rotation speed and strainrate_threshold
samplefile_name=['C0' num2str(filenum(fi,1))];
omegamatrix=repmat(exprec_num(:,1),1,size(exprec_raw,2));
watercontentmatrix=repmat(exprec_num(1,:),size(exprec_raw,1),1);
% omega_sample=omegamatrix(strcmp(samplefile_name,exprec_raw));
omega_sample=60;
watercontent_sample=watercontentmatrix(strcmp(samplefile_name,exprec_raw));
strainrate_coef=110;
strainrate_thresh=omega_sample/180*pi*0.1*strainrate_coef;
% strainrate_thresh=1;

%% Obtain strain rate data by surface line
load([filedir 'C0' num2str(filenum(fi,1)) '\PIVlab_result_comdiy.mat']);
xreal=xoverall*pixeltorealdimension-CCreal(1);
yreal=(yoverall-1080)*-1*pixeltorealdimension-CCreal(2);
particle_disttoCC=sqrt(xreal.^2+yreal.^2);
strainrate_real=strainrate/HFR*0.5;
strainrate_real(particle_disttoCC>0.09)=nan;
strainrate_real_seg=strainrate_real;

% %-----start----extract strain rate region below the surface line---------
for stri=1:size(xreal,2)-1
    xrange=[xreal(1,stri) xreal(1,stri+1)];
    xrange_cell{stri,1}=[xreal(1,stri) xreal(1,stri+1)];
    pixi=find(xrange(1,1)<surface_x_real&surface_x_real<xrange(1,2));%points_insurfaceline_x_index
    piy=surface_y_real(pixi,1);%points_insurfaceline_y
    piy_cell{stri,1}=surface_y_real(pixi,1);%points_insurfaceline_y
    y_segbound(stri,1)=mean(piy);%use this value to segment the column in y overall
    ybsi=find(yreal(:,stri)>y_segbound(stri,1));%y_below_segbound_index
    ybsi_cell{stri,1}=ybsi;
    if isempty(piy)
        strainrate_real_seg(:,stri)=nan;
    else
        strainrate_real_seg(ybsi,stri)=nan;
    end
end
% %----end-----extract strain rate region below the surface line---------

[strr_bound_tmp]=contours(xreal,yreal,strainrate_real_seg,[strainrate_thresh strainrate_thresh]);
strr_bound_tmp=strr_bound_tmp';
circle_index=find(strr_bound_tmp(:,2)==max(strr_bound_tmp(:,2)));
strr_bound_tmp1=strr_bound_tmp(circle_index+1:circle_index+max(strr_bound_tmp(:,2)),:);

% figure(1);hold on
% [c2,h2]=contourf(xreal,yreal,strainrate_real_seg,[strainrate_thresh(fi,1) strainrate_thresh(fi,1)]);
% set(h2, 'edgecolor','none');
% axis equal;
% figure(2);hold on
% histogram(strainrate_real_seg,'Normalization','Probability');

%% Extract shear region lower boundary based on strain rate: clusterdata function
xws=2*dp;%x window size
steps=1*dp;%moving step size
[strr_botbound]=strainrate_bottombound(strr_bound_tmp1,xws,steps);
flowboundary_in_pixel=[surface_x_real,surface_y_real;strr_botbound]/pixeltorealdimension;
flowboundary_in_pixel(:,2)=(flowboundary_in_pixel(:,2)-1080)*-1;
% %---------extract strain rate bottom boundary based on distance to center-----------
disttoCC_strrbotbound=sqrt(strr_botbound(:,1).^2+strr_botbound(:,2).^2);
strrboundinsideCircle_index=find(disttoCC_strrbotbound<0.07);%0.1 is the real container radius
strr_botbound_incircle=strr_botbound(strrboundinsideCircle_index,:);%surface in circle
% datafilename=[char(expnum) 'strinrate_boundary.mat'];
% save(datafilename,'boundary_in_pixel');

%% fitting and find the strainrate depth
marker_color_type='ko';
xwsfit=6*dp;%x window size
stepsfit=1*dp;%moving step size
[strr_depth,strr_depth_normbydp,maxdepth_norm,swfit,maxdepth_topbotpoints]=strr_depth_fun ...
(strr_botbound_incircle,sicx,sicy,dp,xwsfit,stepsfit);
% datafilename=([char(expnum) '_strainrate_depth.mat']);
% save(datafilename,'strr_depth_normbydp');
% save(datafilename,'maxdepth_norm','-append');
% save(datafilename,'swfit','-append');
%% save data

datafilename=[filedir 'C0' num2str(filenum(fi)) '\C0' num2str(filenum(fi)) '_strr-coe' num2str(strainrate_coef) '_strainrate_depth_2.mat'];
save(datafilename,'strr_depth');
save(datafilename,'strr_depth_normbydp','-append');
save(datafilename,'maxdepth_norm','-append');
save(datafilename,'swfit','-append');
save(datafilename,'flowboundary_in_pixel','-append');
save(datafilename,'strainrate_coef','-append');
save(datafilename,'strainrate_thresh','-append');
save(datafilename,'strr_botbound','-append');
save(datafilename,'strainrate_real_seg','-append');
save(datafilename,'surface_x_real','-append');
save(datafilename,'surface_y_real','-append');
save(datafilename,'xreal','-append');
save(datafilename,'yreal','-append');
save(datafilename,'pixeltorealdimension','-append');
save(datafilename,'CCreal','-append');
save(datafilename,'maxdepth_topbotpoints','-append');

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

function [strr_depth,strr_depth_normbydp,maxdepth_norm,swfit,maxdepth_topbotpoints]=strr_depth_fun ...
(strr_botbound,surface_x,surface_y,dp,xwsfit,stepsfit)
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
    midbot(fiti,1:2)=[swfit(fiti) mean(points_tmp(:,2))];
    %% left and right point based on bottom mid point
    botleft_coor=[midbot(fiti,1)-xwsfit/2 midbot(fiti,2)-xwsfit/2*fitangle_a(fiti,1)];
    botright_coor=[midbot(fiti,1)+xwsfit/2 midbot(fiti,2)+xwsfit/2*fitangle_a(fiti,1)];
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
    midtop(fiti,1:2)=[mean(betweenlines_pointscoor(:,1)) mean(betweenlines_pointscoor(:,2))];
    %% flow region depth; dist between midbot and midtop
    if (fitangle_a(fiti,1)+pi/2)<=pi && (fitangle_a(fiti,1)+pi/2)>=pi/2
        strr_depth(fiti,1)=sqrt(sum((midtop(fiti,:)-midbot(fiti,:)).^2));
    else
        strr_depth(fiti,1)=0;
    end
        
    end
strr_depth_normbydp=strr_depth/dp;
maxdepth_norm=max(strr_depth_normbydp);
maxdepth_norm_index=find(strr_depth_normbydp==maxdepth_norm);
maxdepth_topbotpoints=[midtop(maxdepth_norm_index,:);midbot(maxdepth_norm_index,:)];

% figure(1);hold on;plot(swfit,fitangle_a/pi*180,marker_color_type);
% figure(3);hold on;plot(swfit,strr_depth_normbydp,marker_color_type);
end
