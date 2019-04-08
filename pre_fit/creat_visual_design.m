close all
clear all
v = VideoReader('/Users/martin/Desktop/dm.mp4');
stim = uint8(zeros(281,500,120));
ini_ratio = 5/100;
size_fix = 20;
for t = 1:119
    im = read(v,t);
    
    % take out fixation target
    if (t >= 22 && t <= 24) || (t >= 41 && t <= 43) || (t >= 78 && t <= 80) || (t >= 96 && t <= 98)
        val_fix = 255;
    else
        val_fix = 0;
    end
	im(1080/2-size_fix:1080/2+size_fix,1920/2-size_fix:1920/2+size_fix,:) = uint8(val_fix);
    im = imresize(im,ini_ratio,'bilinear','Antialiasing',0);
    
    im(im >= 191 & im <= 191) = uint8(0);
    im(im ~= 0) = uint8(255);
    im(im ~= 0 & im ~= 255) = uint8(0);
    
    im = imresize(im,[281,500],'bilinear','Antialiasing',0);
    stim(:,:,t) = im(:,:,1);
end
stim(:,:,120) = stim(:,:,119);
stim = single(stim);
save('/Users/martin/disks/ae_S/2018/visual/nprf_indiv/raw_data/vis_design.mat','stim');

