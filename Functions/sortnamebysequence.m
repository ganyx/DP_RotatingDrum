function [d]=sortnamebysequence(d)
nameCell = cell(size(d,1),1);
for nc = 1:size(d,1)
%     disp(d(nc).name);
    nameCell{nc} = d(nc).name;
end
addpath('C:\Users\jeffr\Dropbox (Sydney Uni)\projects\2d_packing\data_analysis_code');
d1 = natsort(nameCell);
d = cell2struct(d1,'name',2);