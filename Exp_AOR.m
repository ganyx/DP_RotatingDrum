clc;
clear;
%%
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis')
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\function')
read_filedir='D:\EXP RoDrtest\Exp\for test\original_photos';
save_filedir='D:\EXP RoDrtest\Exp\for test\sharpen_segmented_photos';
figureinfo=dir([read_filedir '\*.jpg']);
% %----------------- C0055_CR is the data: coordinates of rotation drum center (CC), radius of rotation drum (CR)
CCR=load('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\C0055_CR.mat');
CC=CCR.CC;
CC(1,2)=(CC(1,2)-1080)*-1;
CR=CCR.calied_pixelR;
% &----------------- end explanation C0055_CR

for figi=1
%     :size(figureinfo,1)
datafile=figureinfo(figi,1).name;
Expfig=imread([read_filedir '\' datafile]);

%% segment image to only show particle region (this still involve some nonparticle parts)
[BW,segmentedImage] = Exp_segmentImage_fun(Expfig);
BW=double(BW);
% figure(1)
% imshow(BW);
% figure(2)
% imshow(segmentedImage);

%% Obtain AOR
% % -----Obtain pixels in the particle region in using segmentedBW figure
yrep=(repmat([1:1080]',1,1920)-1080)*-1;
xrep=repmat([1:1920],1080,1);
pixdist_toCC=sqrt((xrep-CC(1,1)).^2+(yrep-CC(1,2)).^2); %pixel distance to CC
% figure(3)
% contour(xrep,yrep,pixdist_toCC)
% axis equal
figure(2)
contour(xrep,yrep,BW);hold on;
% imshow(BW);
% set(gca, 'ydir', 'reverse');
BW(pixdist_toCC>CR)=0; % BW: remove pixels outside the circle
BWxcol=xrep(BW==1);
BWycol=yrep(BW==1);
BW_toCCdist=pixdist_toCC(BW==1);
% % ----------boundary method to extract upper surface pixels-----------
boundk=boundary(BWxcol,BWycol,0.9);
% plot(BWxcol,BWycol,'b*',BWxcol(boundk),BWycol(boundk),'r-','LineWidth',2);
% set(gca,'Ydir','reverse');
pixelx_bound=BWxcol(boundk);%boundary pixels
pixely_bound=BWycol(boundk);
boundtoCCdist=BW_toCCdist(boundk);%boundary pixels to CC distance
upsurx=pixelx_bound(boundtoCCdist<CR*0.99);%upper surface particles
upsury=pixely_bound(boundtoCCdist<CR*0.99);

        %% ---------------fitting DAOR and locate sampling area-------------------
        % Set up fittype and options.
        ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.395894484205228 0.936450101590675];
        fitangle_a=[];
        fitangle_b=[];
        num_points_to_fit=100;
        fitdatamidx=[];
        fitdatamidy=[];
        for fiti=1:size(upsurx,1)-num_points_to_fit+1
        [fitresult, gof] = fit(upsurx(fiti:fiti+num_points_to_fit-1),upsury(fiti:fiti+num_points_to_fit-1),ft,opts);
        fitdatamidx(fiti,1)=upsurx(fiti-1+num_points_to_fit/2);
        fitdatamidy(fiti,1)=upsury(fiti-1+num_points_to_fit/2);
        fitangle_a(fiti,1)=atan(fitresult.a);%dynamic angle of repose in degree
        fitangle_b(fiti,1)=fitresult.b;
        gof_fitangle(fiti,1)=gof.rsquare;
        end

        % % -------plot fiting line ------------------
        figure(2)
        hold on
        xft=[0:1:1920];
        yft=tan(max(fitangle_a))*xft+fitangle_b(fitangle_a==max(fitangle_a));
        plot(xft,yft,'r-','LineWidth',2);
        axis equal;
        xlim([0 1920]);
        ylim([0 1080]);
        figure(2);
        hold on;
        %% --------plot max angle point-------------
        coord_maxanglex=fitdatamidx(fitangle_a==max(fitangle_a));
        coord_maxangley=fitdatamidy(fitangle_a==max(fitangle_a));
        plot(coord_maxanglex,coord_maxangley,'go','MarkerSize',10,'LineWidth',1);
        %% --------plot AOR fitting angle transition-------
        figure(3)
        angle_a=fitangle_a/pi*180;
        plot(fitdatamidx,angle_a);
        hold on
        plot(fitdatamidx(angle_a==max(angle_a)),max(angle_a),'go','MarkerSize',10);
        ylabel('Angle(^o)');
        %% ----------plot perpendicular line of the fitting line
        currslope=tan(max(fitangle_a));
        perpslope=-1/currslope;
        perp_b=fitdatamidy(fitangle_a==max(fitangle_a))-(perpslope*fitdatamidx(fitangle_a==max(fitangle_a)));
        figure(2); hold on;
        yft_perp=perpslope*xft+perp_b;
        plot(xft,yft_perp,'g-','LineWidth',1);
        %% ----------plot intersection between perpendicular line and circle
        figure(2);hold on;
        [interx,intery] = linecirc(perpslope,perp_b,CC(1,1),CC(1,2),CR);
        plot(interx(1),intery(1),'y*','MarkerSize',10,'LineWidth',1);
        dist_OI=sqrt((coord_maxanglex-interx(1)).^2+(coord_maxangley-intery(1)).^2);
        %% ----------plot parallel line of fitting line cross intersection-----
        para_b=intery(1)-(currslope*interx(1));
        figure(2); hold on;
        yft_para=currslope*xft+para_b;
        plot(xft,yft_para,'r--','LineWidth',1);
        %% ----------plot rotated figure---------
        % % the rotation is based on the max angle point
        % % the rotation angle is the slope of the max slope
        rotation_angle=(-max(fitangle_a)/pi*180);%minus is anticlockwise
        BWrotate=rotateAround(BW, coord_maxangley*-1+1080, coord_maxanglex,rotation_angle);
        figure(6);imshow(BWrotate);
        hold on
        plot(fitdatamidx(fitangle_a==max(fitangle_a)),(fitdatamidy(fitangle_a==max(fitangle_a))*-1)+1080,'go','MarkerSize',10,'LineWidth',1);
        %% ---------plot sampling region---------
        figure(6);hold on;
        pix_one=6*5;%the assumed pixel for one particle
        S1(1,1)=coord_maxanglex+pix_one;S1(1,2)=coord_maxangley;
        S2(1,1)=coord_maxanglex-pix_one;S2(1,2)=coord_maxangley;
        S3(1,1)=coord_maxanglex+pix_one;S3(1,2)=coord_maxangley-dist_OI;
        S4(1,1)=coord_maxanglex-pix_one;S4(1,2)=coord_maxangley-dist_OI;
        S_coord=[S1;S2;S4;S3];
        plot(S_coord(:,1),S_coord(:,2)*-1+1080,'o','MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
        plot(polyshape(S_coord(:,1),S_coord(:,2)*-1+1080));
        %% obtain velocity field data from sampling region
%         load('E:\EXP RoDrtest\Exp\for test');

%% sharpen image
% % image = gpuArray(imread('E:\EXP RoDrtest\Exp\for test\C0055.MP4_20220225_094202.186.jpg'));
% dimage = im2double(Expfig);
% gradient = convn(dimage,ones(3)./9,'same') - convn(dimage,ones(5)./25,'same');
% amount = 5;
% sharpened = dimage + amount.*gradient;
% % imshow(imresize([dimage, sharpened],0.7));
%% save figure
% cd(save_filedir);
% imwrite(sharpened,[num2str(figi) '.jpg']);
end

