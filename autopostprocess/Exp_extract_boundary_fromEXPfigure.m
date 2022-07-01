clc;
clear;
% % batch extract boundary of particle zone in experiments
%% read file dir
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\github_repo');%for export figure

filenum=[149:155,157:175]';
for fi=1:size(filenum,1);
    tic
filedir=(['D:\EXP RoDrtest\Exp\sharpened\C0' num2str(filenum(fi)) '\']);
iminfo=dir([filedir '*.jpg']);
[iminfo]=sortnamebysequence(iminfo);
[Rtxt,Ctxt]=size(iminfo);
%% locate drum centre and calibration parameters
maskdir=(['C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\pivlab_masks\']);
% load('E:\EXP RoDrtest\Exp\PIVlab_mask.mat');
if filenum(fi)>156
    load([maskdir 'PIVlab_mask1_c0157-c0162.mat'])
else
    load([maskdir 'PIVlab_mask1_c0149-c0155.mat'])
end

x=cell2mat(maskiererx(1,1));
y=cell2mat(maskierery(1,1));
% drum_cen=[mean(x) mean(y)];

p1=[x(1) y(1)];
p2=[x(2) y(2)];
p3=[x(3) y(3)];

%% calibration para
figure(2)
fig=imread([filedir 'C0' num2str(filenum(fi)) '_00001.jpg']);
imshow(fig)
if filenum(fi)>156
    pixeltorealdimension=0.10955/835;%c0157-c0162
else
    pixeltorealdimension=0.10955/820;%c0149-c0155
end

% % ------------------drum center------------------
hold on
[CC,Radius]=CircleThru3Dots(p1,p2,p3);
% % --------------------calibration coefficient-----------------
% calicoe=0.96; %calibration coefficient

theta=linspace(0,2*pi,1001);

%% plot calibration
% hold on;
% xcorr1=CC(1)+0.1/cali*cos(theta);
% ycorr1=CC(2)+0.1/cali*sin(theta);
% plot(xcorr1,ycorr1,'r--');

xcorr2=CC(1)+0.1/pixeltorealdimension*cos(theta);
ycorr2=CC(2)+0.1/pixeltorealdimension*sin(theta);
plot(xcorr2,ycorr2,'r--');

xlim([0 1920]);ylim([0 1080]);
% close figure 2
%% locate boundary
bx=[];
by=[];
processquant=100;
for i=1:processquant
    process=i/processquant
        image=imread([filedir iminfo(i).name]);
    % % 1. run Exp_strainrate_pivlab
    % % 2. load the image manuely e.g. C0083_00100
    [BW,maskedImage] = Exp_segmentImage_fun_grey_c0154(image);
    xtmp=[1:1920];
    ytmp=flip([1:1080]');
    xtmprepmat=repmat(xtmp,1080,1);
    ytmprepmat=repmat(ytmp,1,1920);
    double(BW);
    ztmpreshape=reshape(double(BW),[],1);
    % figure(4)
    % contourf(xtmprepmat,ytmprepmat,double(BW));
    xtmpreshape=reshape(xtmprepmat,[],1);
    ytmpreshape=reshape(ytmprepmat,[],1);
    CCreal=[CC(1)*pixeltorealdimension,((CC(2)-1080)*pixeltorealdimension)*-1];
    disttoCCreal=sqrt((xtmprepmat*pixeltorealdimension-CCreal(1)).^2+(ytmprepmat*pixeltorealdimension-CCreal(2)).^2);
    disttoCCrealreshape=reshape(disttoCCreal,[],1);
    index_incircle=find(disttoCCrealreshape<=0.1);
    index_particle=find(ztmpreshape==1);
    index_particle_incircle=index_incircle(ismember(index_incircle,index_particle));
    x_particle_incircle=xtmpreshape(index_particle_incircle)*pixeltorealdimension;
    y_particle_incircle=ytmpreshape(index_particle_incircle)*pixeltorealdimension;
    % figure(5)
    % plot(x_particle_incircle,y_particle_incircle,'ko');
    % hold on;
    k=boundary(x_particle_incircle,y_particle_incircle,0.8);
    tmpx=x_particle_incircle(k);
    tmpy=y_particle_incircle(k);
    tmpdist=sqrt((tmpx-CCreal(1)).^2+(tmpy-CCreal(2)).^2);
    index_obtain_upper_bound=find(tmpdist<0.09);
    figure(2)
    hold on
    xplot=tmpx(index_obtain_upper_bound)/pixeltorealdimension;
    yplot=tmpy(index_obtain_upper_bound)/pixeltorealdimension;
    plot(xplot,(yplot-1080)*-1,'ro');
    hold on
%     xlim([0.05 0.3]);
    axis equal
    ylim([0 1080]);
    bx{i,1}=tmpx(index_obtain_upper_bound);
    by{i,1}=tmpy(index_obtain_upper_bound);

end
datafilename=[filedir 'C0' num2str(filenum(fi)) '_particleboundary_2.mat'];
save(datafilename,'bx');
save(datafilename,'by','-append');
cd(filedir);
export_fig uppersurface_cloudpoints_2.jpg -native
close figure 2
toc
end