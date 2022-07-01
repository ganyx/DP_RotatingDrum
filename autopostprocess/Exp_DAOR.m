clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('D:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
filenum=[149:155,157:172]';
%%
for fi=1:size(filenum,1)
    disp(['C0' num2str(filenum(fi))]);
% figureinfo=dir([read_filedir '\*.jpg']);
%% extract calibration and Container center info
if filenum<156
    CCR=load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\C0149-C0155_CR.mat');
    load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\calibration_coe_c0149-c0155.mat');
else
    CCR=load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\C0157-C0172_CR.mat');
    load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\calibration_coe_c0157-c0172.mat');
end
CCreal=CCR.CC*pixeltorealdimension;
CCreal(1,2)=(CCreal(1,2)-1080*pixeltorealdimension)*-1;

%% load upper surface data surface_x surface_y
load([filedir 'C0' num2str(filenum(fi,1)) '\C0' num2str(filenum(fi,1)) '_surfacepoints_2.mat']);

%% extrac surface points inside container
disttoCC=sqrt((surface_x-CCreal(1)).^2+(surface_y-CCreal(2)).^2);
surfinsideCircle_index=find(disttoCC<0.07);%0.1 is the real container radius
sicx=surface_x(surfinsideCircle_index);%surface in circle
sicy=surface_y(surfinsideCircle_index);

%% Obtain DAOR
    %% ---------------fitting DAOR and locate sampling area-------------------
    % Set up fittype and options.
    ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.395894484205228 0.936450101590675];
    fitangle_a=[];
    fitangle_b=[];
    num_points_to_fit=200;
    interval=floor(size(sicx,1)/num_points_to_fit);
    intervalpara=5;
    samples=[1:num_points_to_fit/intervalpara:size(sicx,1)]';
    fitdatamidx=[];
    fitdatamidy=[];
    gof_fitangle=[];
    for fiti=1:size(samples,1)-intervalpara
    [fitresult, gof] = fit(sicx(samples(fiti):samples(fiti)+num_points_to_fit),sicy(samples(fiti):samples(fiti)+num_points_to_fit),ft,opts);
    fitdatamidx(fiti,1)=mean(sicx(samples(fiti):samples(fiti)+num_points_to_fit));
    fitdatamidy(fiti,1)=mean(sicy(samples(fiti):samples(fiti)+num_points_to_fit));
    fitangle_a(fiti,1)=atan(fitresult.a);%dynamic angle of repose in degree
    fitangle_b(fiti,1)=fitresult.b;
    gof_fitangle(fiti,1)=gof.rsquare;
    end

% % -------plot fiting line ------------------
DAOR=max(fitangle_a);
% figure(2)
% hold on
% xft=[0:1:1920]*pixeltorealdimension;
% yft=tan(max(fitangle_a))*xft+fitangle_b(fitangle_a==max(fitangle_a));
% plot(xft,yft,'r-','LineWidth',2);
% axis equal;
% % xlim([0 1920]);
% % ylim([0 1080]);
% % figure(2);
% % hold on;
%% --------plot max angle point-------------
% coord_maxanglex=fitdatamidx(fitangle_a==max(fitangle_a));
% coord_maxangley=fitdatamidy(fitangle_a==max(fitangle_a));
% plot(coord_maxanglex,coord_maxangley,'go','MarkerSize',10,'LineWidth',1);
%% --------plot AOR fitting angle transition-------
% figure(3)
% angle_a=fitangle_a/pi*180;
% plot(fitdatamidx,angle_a);
% hold on
% plot(fitdatamidx(angle_a==max(angle_a)),max(angle_a),'go','MarkerSize',10);
% ylabel('Angle(^o)');
%% save data

datafilename=[filedir 'C0' num2str(filenum(fi)) '\C0' num2str(filenum(fi)) '_DAOR_2.mat'];
save(datafilename,'fitangle_a');
save(datafilename,'fitangle_b','-append');
save(datafilename,'gof_fitangle','-append');
save(datafilename,'fitdatamidx','-append');
save(datafilename,'fitdatamidy','-append');
save(datafilename,'DAOR','-append');
save(datafilename,'num_points_to_fit','-append');
end


