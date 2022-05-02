clc;
clear;
% plot strain rate from results of PIVlab
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
filedir=(['e:\EXP RoDrtest\Exp\analysis data\']);
fileinfo=dir([filedir '*.mat']);
[fileinfo]=sortnamebysequence(fileinfo);
[Rtxt,Ctxt]=size(fileinfo);

for fi=5
%     :Rtxt
    %% load data
    data = load([filedir fileinfo(fi).name]);
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
                tmpv=cell2mat(cellv(si));
                tmpstrr=cell2mat(cellstrr(si));
                tmpsher=cell2mat(cellsher(si));
%             else
                utotal=utotal+tmpu;
                vtotal=vtotal+tmpv;
                strrtotal=strrtotal+tmpstrr;
                shertotal=shertotal+tmpsher;
%             end
        end

     uave=utotal/totalframes/1000;
     vave=vtotal/totalframes/1000;
%       strrave=strrtotal/totalframes*1000000;
     strrave1=strrtotal/totalframes/1000;
     sherave1=shertotal/totalframes/1000;
    %% ========= plot strain rate contour: accumulate after strain rate each frame ===========
    figure(1)
%     strrave(strrave>10000)=10000;
%     strrave(strrave<-10000)=-10000;
    [c1,h1]=contourf(cell2mat(cellx(si)),(cell2mat(celly(si))-max(cell2mat(celly(si)),[],'all'))*-1,sherave1,200);
    set(h1, 'edgecolor','none');
%     set(h,'Ydir','reverse');
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