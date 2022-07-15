clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
filedir=('D:\EXP RoDrtest\Exp\sharpened\');
dp=0.002;%particle diameter
filenum=[149:155,157:174]';
[exprec_num,exprec_txt,exprec_raw]=xlsread('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\parametric_study.xlsx',8);
omega=exprec_num(2:7,1);
watercontent=exprec_num(1,2:5);
%%
flowdepthmatrix=zeros(size(omega,1),size(watercontent,2));
% color=colormap(turbo(16));
color=flip(colormap(hot(16)));
for fi=2:size(exprec_raw,1)
    for fj=2:size(exprec_raw,2)
        load([filedir exprec_raw{fi,fj} '\' exprec_raw{fi,fj} '_strr-coe50_strainrate_depth_2.mat']);
        flowdepthmatrix(fi-1,fj-1)=maxdepth_norm;
        figure(1);hold on;
%         scatter(omega(fi-1,1)/60,flowdepthmatrix(fi-1,fj-1),'o','filled','MarkerFaceColor',color(fj*2,:));
        plot(omega(fi-1,1)/60,flowdepthmatrix(fi-1,fj-1),'o', ...
            'MarkerFaceColor',color(fj*3,:), ...
            MarkerEdgeColor='none');
    end
end

xlabel('Rotation speed (rpm)');
ylabel('Normalized flow region depth (D/d)');
set(gca,'FontSize',18);
set(gca,'FontName','Times New Roman');
legend('dry', 'm=0.1','m=0.2','m=0.5','Location','northwest');
set(gca,'LineWidth',1);
box on
xlim([0 2.6]);