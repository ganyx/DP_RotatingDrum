I=imread('G:\EXP RoDrtest\Exp\for test\original_photos\C0055.MP4_20220225_094202.186.jpg');
im=rgb2gray(I);
im=imadjust(im);
im=im2bw(im,240/255);
im=medfilt2(im);

 figure,imshow(im);
Rmin = 30; 
Rmax = 80;
im=~im;
[center, radius] = imfindcircles(im,[Rmin Rmax],'Sensitivity',0.9,'ObjectPolarity','bright');
viscircles(center,radius);

cir=size(center);
tot_cir=cir(1);

 im=imopen(im,strel('disk',2));
im=imfill(im,'holes')
% figure,imshow(im);

[Ilabel num] = bwlabel(im);
disp(num);
Iprops = regionprops(Ilabel);
Ibox = [Iprops.BoundingBox]; 
Ibox = reshape(Ibox,[4 num]);


b=bwboundaries(im);

a=size(b);
disp('Number of Bolts=');
disp(a(1)-tot_cir);
bolt=a(1)-tot_cir;
disp('Number of Nuts= ');
disp(tot_cir);

figure,imshow(I);
%figure,imshow(im);

for cnt = 1:num 
    rectangle('position',Ibox(:,cnt),'edgecolor','r');
end
title(sprintf('No of Bolts= %1.0f , No of Nuts= %1.0f',bolt,tot_cir));



 