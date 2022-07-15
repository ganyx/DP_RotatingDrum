clc;
clear;
% % batch extract boundary of particle zone in experiments
%% read file dir
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
filedir=(['d:\EXP RoDrtest\Exp\sharpened\']);
filenum=[149:154,157:174]';%149:154,
%% locate drum centre and calibration parameters
for fi=1:size(filenum,1)
    picinfo=dir([filedir 'C0' num2str(filenum(fi,1)) '\*.jpg']);
    [picinfo]=sortnamebysequence(picinfo);
    [Rtxt,Ctxt]=size(picinfo);
    load([filedir 'C0' num2str(filenum(fi,1)) '\C0' num2str(filenum(fi,1)) '_strr-coe80_strainrate_depth_2.mat']);
    cd('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\Roration drum experiment\Dataanalysis\results\video_surf-strr_bound');    
    video = VideoWriter(['c0' num2str(filenum(fi,1)) '_strr-coe80_surf-strr_bound']); %create the video object 
    open(video); %open the file for writing
        for i=1:20

            image=imread([filedir 'C0' num2str(filenum(fi,1)) '\' picinfo(i).name]);
            figure(1)
            imshow(image);
            hold on;
            boundary_in_pixel_vel=[surface_x_real,surface_y_real;strr_botbound]/pixeltorealdimension;
            %---------plot surface & strainrate boundary line--------------
            plot(boundary_in_pixel_vel(:,1)+(CCreal(1)/pixeltorealdimension),(boundary_in_pixel_vel(:,2)-1080)*-1-(CCreal(2)/pixeltorealdimension),'r.');
            %---------plot max depth points--------------
%             plot((maxdepth_topbotpoints(:,1)+CCreal(1))/pixeltorealdimension,((maxdepth_topbotpoints(:,2)-CCreal(2))/pixeltorealdimension-1080)*-1,'ro');
            plot((maxdepth_topbotpoints(:,1)+CCreal(1))/pixeltorealdimension, ...
                ((maxdepth_topbotpoints(:,2)+CCreal(2))/pixeltorealdimension-1080)*-1,'go-', ...
                'LineWidth',2, ...
                'MarkerSize',6);
            title(['C0' num2str(filenum(fi,1))]);
            %% save as video
            F = getframe(gcf);
            [X, Map] = frame2im(F);
            writeVideo(video,X); %write the image to file
        %% save as gif
%             F=getframe(gcf);
%             I=frame2im(F);
%             [I,map]=rgb2ind(I,256);
%             if i == 1
%                 imwrite(I,map,['C0' num2str(filenum(fi,1)) '_strthre' num2str(round(strainrate_thresh)) '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.2);
%             else
%                 imwrite(I,map,['C0' num2str(filenum(fi,1)) '_strthre' num2str(round(strainrate_thresh)) '.gif'],'gif','WriteMode','append','DelayTime',0.2);
%             end
%      
%             close figure 1
        end
        close(video); %close the file
        clear image video
end
% video = VideoWriter('yourvideo.avi'); %create the video object
% open(video); %open the file for writing
% for ii=1:N %where N is the number of images
%   I = imread('the ith image.jpg'); %read the next image
%   writeVideo(video,I); %write the image to file
% end
% close(video); %close the file