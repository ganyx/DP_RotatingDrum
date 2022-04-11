clc;
clear;
%% read file dir
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\unsaturated\dataanalysis\');
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code\functions\');
filedir=(['E:\EXP RoDrtest\Exp']);
videoinfo=dir([filedir '\*.mp4']);
[videoinfo]=sortnamebysequence(videoinfo);
[Rtxt,Ctxt]=size(videoinfo);
%%
for fi=35:38
%     :Rtxt
    videoname=videoinfo(fi).name;
    vr = VideoReader(['E:\EXP RoDrtest\Exp\' videoname]);
    videoNo=split(videoname,'.');
    numFrames = ceil(vr.FrameRate*vr.Duration);
    cd(['E:\EXP RoDrtest\Exp\sharpened\' videoNo{1}]);
    for fri=1:numFrames
%         :10
%         :numFrames
        frame = read(vr,fri);
        frame_gray=rgb2gray(frame);
        %% sharpen image method 1

        dimage = im2double(frame_gray);
        gradient = convn(dimage,ones(3)./9,'same') - convn(dimage,ones(5)./25,'same');
        amount = 5;
        sharpened = dimage + amount.*gradient;

%     imshow(imresize([dimage, sharpened],0.7)); %compare original and sharpened image
        %% sharpen image method 2
%         sharpened = imsharpen(frame,'Radius',2,'Amount',5);
%         imshow(imresize([frame, sharpened],0.7)); %compare original and sharpened image

        %% write image        
        imwrite(sharpened,[videoNo{1} '_' num2str(fri,'%05.f') '.jpg'])
    end
    close all;
end