ims = dir('..\Output\');
n = length(ims) - 2;
images = cell(n,1);

% Start at 3 because for some reason ims includes . and ..
for i=1:n
    images{i} = imread(strcat(strcat('..\Output\',sprintf('%d',i)),'.png'));
end

v = VideoWriter('..\results5.mp4','MPEG-4');
open(v);

for i=1:n
    writeVideo(v,images{i});
end
close(v);