clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('d:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
filenum=[149:155,157:174]';
[exprec_num,exprec_txt,exprec_raw]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',8);
omega=exprec_num(2:7,1);
watercontent=exprec_num(1,2:5);
%%
DAORmatrix=zeros(size(omega,1),size(watercontent,2));
color=flip(colormap(hot(16)));
for fi=2:size(exprec_raw,1)
    for fj=2:size(exprec_raw,2)
        load([filedir exprec_raw{fi,fj} '\' exprec_raw{fi,fj} '_DAOR_2.mat']);
        DAORmatrix(fi-1,fj-1)=DAOR/pi*180;
        figure(1);hold on;
        scatter(omega(fi-1,1)/60,DAORmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*3,:));
    end
end

xlabel('Rotation speed (rpm)');
ylabel('Dynamic angle of repose (^\circ)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('dry', 'm=0.1','m=0.2','m=0.5');
set(gca,'LineWidth',1);
box on
% xlim([0 160]);