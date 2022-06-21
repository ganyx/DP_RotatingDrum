clc;
clear;
% % batch extract boundary of particle zone in experiments
%% read file dir
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function');
filedir=(['E:\EXP RoDrtest\Exp\sharpened\C0084\']);

videoinfo=dir([filedir '*.jpg']);
[videoinfo]=sortnamebysequence(videoinfo);
[Rtxt,Ctxt]=size(videoinfo);
%% locate drum centre and calibration parameters
load('E:\EXP RoDrtest\Exp\PIVlab_mask.mat');
x=cell2mat(maskiererx(1,1));
y=cell2mat(maskierery(1,1));

% drum_cen=[mean(x) mean(y)];

p1=[x(1) y(1)];
p2=[x(2) y(2)];
p3=[x(3) y(3)];
cali=0.236/1295.64;
pixeltorealdimension=cali/0.96;
% % ------------------drum center------------------
hold on
[CC,Radius]=CircleThru3Dots(p1,p2,p3);
% % --------------------calibration coefficient-----------------
calicoe=0.96; %calibration coefficient
% % figure(1)
% theta=linspace(0,2*pi,1001);
% hold on;
% xcorr1=CC(1)+0.1/cali*cos(theta);
% ycorr1=CC(2)+0.1/cali*sin(theta);
% plot(xcorr1,ycorr1,'r--');
% 
% xcorr2=CC(1)+0.1/cali*calicoe*cos(theta);
% ycorr2=CC(2)+0.1/cali*calicoe*sin(theta);
% plot(xcorr2,ycorr2,'y--');
%% locate boundary
bx=[];
by=[];
processquant=100;
for i=1:processquant
process=i/processquant
    image=imread([filedir videoinfo(i).name]);
% % 1. run Exp_strainrate_pivlab
% % 2. load the image manuely e.g. C0083_00100
[BW,maskedImage] = Exp_segmentImage_fun_grey_c0084(image);
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
CCreal=CC*pixeltorealdimension;
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
index_obtain_upper_bound=find(tmpdist<0.095);
plot(tmpx(index_obtain_upper_bound),tmpy(index_obtain_upper_bound),'r-');
hold on
xlim([0.05 0.3]);
axis equal
ylim([0 0.2]);
bx{i,1}=tmpx(index_obtain_upper_bound);
by{i,1}=tmpy(index_obtain_upper_bound);

end
datafilename=['C0084' '_particleboundary.mat'];
save(datafilename,'bx');
save(datafilename,'by','-append');