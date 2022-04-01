clc;
clear;

fig=imread('E:\EXP RoDrtest\Exp\C0055.MP4_20220223_205621.535.jpg');
figure(1)
imshow(fig);

load('E:\EXP RoDrtest\Exp\PIVlab_mask.mat');
x=cell2mat(maskiererx(1,1));
y=cell2mat(maskierery(1,1));

% drum_cen=[mean(x) mean(y)];

p1=[x(1) y(1)];
p2=[x(2) y(2)];
p3=[x(3) y(3)];
cali=0.236/1295.64;

%%
hold on
[CC,Radius]=CircleThru3Dots(p1,p2,p3);

%%
calicoe=0.96; %calibration coefficient
figure(1)
theta=linspace(0,2*pi,1001);
hold on;
xcorr1=CC(1)+0.1/cali*cos(theta);
ycorr1=CC(2)+0.1/cali*sin(theta);
plot(xcorr1,ycorr1,'r--');

xcorr2=CC(1)+0.1/cali*calicoe*cos(theta);
ycorr2=CC(2)+0.1/cali*calicoe*sin(theta);
plot(xcorr2,ycorr2,'y--');

%%
% % ----------make the cell of mask for PIVlab readable format-----------
markouter=[0,0;0,1080;1920,1080;1920,0;0,0];
maskiererx{1,1}=[markouter(:,1);xcorr2'];
maskiererx{1,2}=[markouter(:,1);xcorr2'];
maskierery{1,1}=[markouter(:,2);ycorr2'];
maskierery{1,2}=[markouter(:,2);ycorr2'];
calied_pixelR=0.1/cali*calicoe;
% % ----------save mask and circle data
% maskfilename=['C0055_mask.mat'];
% save(maskfilename,'maskiererx');
% save(maskfilename,'maskierery','-append');
centrefilename=['C0055_CR.mat'];
save(centrefilename,'CC');
save(centrefilename,'calied_pixelR','-append');

%%
function [CC,Radius]=CircleThru3Dots(A,B,C)
Ah=A*A';
Bh=B*B';
Ch=C*C';
CC=zeros(size(A));
G=(C(2)-B(2))*A(1)+(A(2)-C(2))*B(1)+(B(2)-A(2))*C(1);
CC(1)=((Bh-Ch)*A(2)+(Ch-Ah)*B(2)+(Ah-Bh)*C(2))/(2*G);
CC(2)=-((Bh-Ch)*A(1)+(Ch-Ah)*B(1)+(Ah-Bh)*C(1))/(2*G);
Radius=sqrt((A-CC)*(A-CC)');
theta=linspace(0,2*pi,101);
x=CC(1)+Radius*cos(theta);
y=CC(2)+Radius*sin(theta);
plot(x,y,'r-')
ABC=[A;B;C];
hold on
plot(ABC(:,1),ABC(:,2),'b.','markersize',20)
plot(CC(1),CC(2),'r.','markersize',20)
grid on
box off
axis equal
end

