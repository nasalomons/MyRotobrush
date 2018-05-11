
ims = dir('..\Input\');
n = length(ims);
images = cell(n,1);

% Start at 3 because for some reason ims includes . and ..
for i=3:n
    images{i} = imread(fullfile('..\Input',ims(i).name));
end
images = images(3:n,1);
n = n - 2;
image = imread('..\Input\1.jpg');

[polymask, xi, yi] = roipoly(image);
outline = bwmorph(polymask,'remove');
[r,c] = find(outline);

windowindeces = cell(round(size(r,1)/20)-1,3);
[cx, cy, cs] = improfile(rgb2gray(image),xi,yi,round(size(r,1)/20));
for i=1:(size(cs,1)-1)
    if round(cy(i,1)) <= 25 
        cy(i,1) = 26;
    end
    if round(cy(i,1)) >= size(image,1) - 25
        cy(i,1) = size(image,1) - 26;
    end
    if round(cx(i,1)) <= 25
        cx(i,1) = 26;
    end
    if round(cx(i,1)) >= size(image,2) - 25
        cx(i,1) = size(image,2) - 26;
    end
    windowindeces{i,3} = [round(cy(i,1)) round(cx(i,1))];
    windowindeces{i,1} = [round(cy(i,1))-25 round(cx(i,1))-25];
    windowindeces{i,2} = [round(cy(i,1))+25 round(cx(i,1))+25];
end

windows = cell(size(windowindeces,1),1);
boxes = cell(size(windowindeces,1),1);
for i=1:size(windowindeces,1)
    b = zeros(size(images{1},1),size(images{1},2));
    c = 1;
    m = zeros(51,51,3);
    for j=windowindeces{i,1}(1,1):windowindeces{i,2}(1,1)
        r = 1;
        for k=windowindeces{i,1}(1,2):windowindeces{i,2}(1,2)
            m(r,c,:) = images{1}(j,k,:);
            b(j,k) = 1;
            r = r + 1;
        end
        c = c + 1;
    end
    windows{i,1} = m;
    boxes{i,1} = imbinarize(b);
end

%{
imshow(images{1});
hold on;
for i=1:size(boxes,1)
    curr = boxes{i,1};
    [L,temp] = bwlabel(curr);
    box = regionprops(L, 'BoundingBox');
    rectangle('Position',box(1).BoundingBox,'EdgeColor','yellow');
end
hold off;
saveas(temp,strcat('..\',sprintf('%d',1)),'png');


temp = imshow(image); axis on;
hold on;
[L,~] = bwlabel(boxes{1});
box = regionprops(L, 'BoundingBox');
rectangle('Position',box(1).BoundingBox,'EdgeColor','yellow');
hold off;
saveas(temp,strcat('..\',sprintf('%d',2)),'png');
%}

%shrink/thicken to make buffer b/e fg/bg
foreground = bwmorph(polymask,'shrink',3);
background = ~(bwmorph(polymask,'thicken',3));
disttrans = bwdist(bwmorph(foreground,'remove'));

windowcms = cell(size(windowindeces,1),1);
oldGFB = cell(size(windowindeces,1),2);
for i=1:size(windowindeces,1)
    localfg = zeros(51,51);
    localbg = zeros(51,51);
    for r=1:size(localfg,1)
        for c=1:size(localfg,2)
            localfg(r,c) = foreground(r-1+windowindeces{i,1}(1,1),c-1+windowindeces{i,1}(1,2));
            localbg(r,c) = background(r-1+windowindeces{i,1}(1,1),c-1+windowindeces{i,1}(1,2));
        end
    end
    
    lab = rgb2lab(image(windowindeces{i,1}(1,1):windowindeces{i,2}(1,1), windowindeces{i,1}(1,2):windowindeces{i,2}(1,2),:));
    l = lab(:,:,1);
    a = lab(:,:,2);
    b = lab(:,:,3);

    gmmf = fitgmdist([double(l(localfg == 1)) double(a(localfg == 1)) double(b(localfg == 1))] ,3,'Options',statset('MaxIter',3000),'RegularizationValue',.000000000001);
    gmmb = fitgmdist([double(l(localbg == 1)) double(a(localbg == 1)) double(b(localbg == 1))] ,3,'Options',statset('MaxIter',3000),'RegularizationValue',.000000000001);

    l = double(reshape(l,[size(l,1)*size(l,2),1]));
    a = double(reshape(a,[size(l,1)*size(l,2),1]));
    b = double(reshape(b,[size(l,1)*size(l,2),1]));

    pdff = pdf(gmmf,[l a b]);
    pdfb = pdf(gmmb,[l a b]);
    
    windowcms{i,1} = pdff ./ (pdff+pdfb);
    windowcms{i,1} = reshape(windowcms{i,1},[size(windows{i,1},1) size(windows{i,1},2)]);
    
    oldGFB{i,1} = gmmf;
    oldGFB{i,2} = gmmb;
end
%{
temp = imshow(windowcms{1});
saveas(temp,strcat('..\',sprintf('%d',8)),'png');
%}


% color confidence
windowcolorconf = cell(size(windows,1),1);
for i=1:size(windowcms,1)
    localfg = foreground(windowindeces{i,1}(1,1):windowindeces{i,2}(1,1), windowindeces{i,1}(1,2):windowindeces{i,2}(1,2));
    localdist = disttrans(windowindeces{i,1}(1,1):windowindeces{i,2}(1,1), windowindeces{i,1}(1,2):windowindeces{i,2}(1,2));

    numerator = 0;
    denominator = 0;
    for r=1:size(windowcms{i,1},1)
        for c=1:size(windowcms{i,1},2)
            L = localfg(r,c);
            p = windowcms{i,1}(r,c);
            d = localdist(r,c);
            w = exp(-(d^2)/(round(size(windowcms{i,1},1)/2)));
            numerator = numerator + (abs(L - p) * w);
            denominator = denominator + w;
        end
    end
    windowcolorconf{i,1} = 1 - (numerator/denominator);
end

% shape models
shapemodels = cell(size(windows,1),1);
for i=1:size(windowindeces,1)
    shapemodels{i,1} = zeros(size(windows{i,1},1), size(windows{i,1},2));
    for j=windowindeces{i,1}(1,1):windowindeces{i,2}(1,1)
        for k=windowindeces{i,1}(1,2):windowindeces{i,2}(1,2)
            L = foreground(j,k);
            d = -(disttrans(j,k))^2;
            a = (size(windowcms{i,1},1) - 2)/((1 - .85)^2);
            if windowcolorconf{i,1} > .85
                s = 2 + (a * (windowcolorconf{i,1} - .85)^2);
            else
                s = 2;
            end
            shapemodels{i,1}(j-windowindeces{i,1}(1,1)+1,k-windowindeces{i,1}(1,2)+1) = 1 - exp(d/(s^2));
        end
    end
end
%{
temp = imshow(foreground(windowindeces{1,1}(1,1):windowindeces{1,2}(1,1), windowindeces{1,1}(1,2):windowindeces{1,2}(1,2)));
saveas(temp,strcat('..\',sprintf('%d',14)),'png');

temp = imshow(shapemodels{1});
saveas(temp,strcat('..\',sprintf('%d',20)),'png');
%}

opticalflow = opticalFlowFarneback;
flow = estimateFlow(opticalflow,rgb2gray(image));

numwindows = size(windows,1);
windowsize = size(windows{1,1},1);

lazybw = lazysnapping(image,superpixels(image,1000),foreground,~foreground);
lazybw = bwmorph(lazybw,'remove');    
temp = image;
R = temp(:,:,1);
G = temp(:,:,2);
B = temp(:,:,3);
R(lazybw) = 255;
G(lazybw) = 0;
B(lazybw) = 0;
temp(:,:,1) = R;
temp(:,:,2) = G;
temp(:,:,3) = B;
temp = imshow(temp);
saveas(temp,strcat('..\Output\',sprintf('%d',1)),'png');
prev = image;

for i=2:size(images)
    prev = image;
    image = imread(strcat(strcat('..\Input\',sprintf('%d',i)),'.jpg'));
    
    [L,~] = bwlabel(foreground);
    box = regionprops(L, 'BoundingBox');
    
    max = 0;
    for j=1:numel(box)
        if box(j).BoundingBox(1,3)*box(j).BoundingBox(1,4) > max 
            rect = box(j).BoundingBox;
            max = box(j).BoundingBox(1,3)*box(j).BoundingBox(1,4);
        end
    end
    
    points1 = detectSURFFeatures(rgb2gray(prev), 'ROI', rect, 'MetricThreshold', 100);
    points2 = detectSURFFeatures(rgb2gray(image), 'ROI', rect, 'MetricThreshold', 100);
    
    [f1, vpoints1] = extractFeatures(rgb2gray(prev),points1);
    [f2, vpoints2] = extractFeatures(rgb2gray(image),points2);
    
    pairs = matchFeatures(f1,f2);
    
    matchedpoints1 = vpoints1(pairs(:,1),:);
    matchedpoints2 = vpoints2(pairs(:,2),:);
    
    tform = estimateGeometricTransform(matchedpoints1,matchedpoints2,'affine');
    
    for j=1:numwindows
        [a, b] = transformPointsForward(tform,windowindeces{j,3}(1,2),windowindeces{j,3}(1,1));
        a = round(a);
        b = round(b);
        if b <= 25 
            b = 26;
        end
        if b >= size(image,1) - 25
            b = size(image,1) - 26;
        end
        if a <= 25
            a = 26;
        end
        if a >= size(image,2) - 25
            a = size(image,2) - 26;
        end
        windowindeces{j,3} = [round(b), round(a)];
        windowindeces{j,1} = [round(b)-25, round(a)-25];
        windowindeces{j,2} = [round(b)+25, round(a)+25];
    end
    
   %{
    figure; showMatchedFeatures(prev,image,matchedpoints1,matchedpoints2,'montage');
    
for l=1:size(windowindeces,1)
    b = zeros(size(images{1},1),size(images{1},2));
    c = 1;
    for j=windowindeces{l,1}(1,1):windowindeces{l,2}(1,1)
        r = 1;
        for k=windowindeces{l,1}(1,2):windowindeces{l,2}(1,2)
            b(j,k) = 1;
            r = r + 1;
        end
        c = c + 1;
    end
    boxes{l,1} = imbinarize(b);
end

temp = imshow(image);
hold on;
for k=1:size(boxes,1)
    curr = boxes{k,1};
    [L,temp] = bwlabel(curr);
    box = regionprops(L, 'BoundingBox');
    rectangle('Position',box(1).BoundingBox,'EdgeColor','yellow');
end
hold off;
saveas(temp,strcat('..\Output\',sprintf('%d',i)),'png');
%}

    flow = estimateFlow(opticalflow,rgb2gray(image));
    X = flow.Vx;
    Y = flow.Vy;
    O = flow.Orientation;
    M = flow.Magnitude;
    
    winflows = cell(numwindows,4);
    for j=1:numwindows
        num = 0;
        winflows{j,1} = 0;
        winflows{j,2} = 0;
        winflows{j,3} = 0;
        winflows{j,4} = 0;
        for r=windowindeces{j,1}(1,1):windowindeces{j,2}(1,1)
            for c=windowindeces{j,1}(1,2):windowindeces{j,2}(1,2)
                if foreground(r,c) == 1      
                    winflows{j,1} = winflows{j,1} + X(r,c);
                    winflows{j,2} = winflows{j,2} + Y(r,c);
                    winflows{j,3} = winflows{j,3} + O(r,c);
                    winflows{j,4} = winflows{j,4} + M(r,c);
                    num = num + 1;
                end
            end
        end
        if num ~= 0
            winflows{j,1} = winflows{j,1} / num;
            winflows{j,2} = winflows{j,2} / num;
            winflows{j,3} = winflows{j,3} / num;
            winflows{j,4} = winflows{j,4} / num;
        end
    end
    
    for j=1:numwindows        
        windowindeces{j,3} = [round(windowindeces{j,3}(1,1)+winflows{j,2}) round(windowindeces{j,3}(1,2)+winflows{j,1})];
        if windowindeces{j,3}(1,1) <= 25 
            windowindeces{j,3}(1,1) = 26;
        end
        if windowindeces{j,3}(1,1) >= size(image,1) - 25
            windowindeces{j,3}(1,1) = size(image,1) - 26;
        end
        if windowindeces{j,3}(1,2) <= 25
            windowindeces{j,3}(1,2) = 26;
        end
        if windowindeces{j,3}(1,2) >= size(image,2) - 25
            windowindeces{j,3}(1,2) = size(image,2) - 26;
        end
        windowindeces{j,1} = [(windowindeces{j,3}(1,1)-25) (windowindeces{j,3}(1,2)-25)];
        windowindeces{j,2} = [(windowindeces{j,3}(1,1)+25) (windowindeces{j,3}(1,2)+25)]; 
    end
    
    
    for l=1:size(windowindeces,1)
        b = zeros(size(images{1},1),size(images{1},2));
        c = 1;
        for j=windowindeces{l,1}(1,1):windowindeces{l,2}(1,1)
            r = 1;
            for k=windowindeces{l,1}(1,2):windowindeces{l,2}(1,2)
                b(j,k) = 1;
                r = r + 1;
            end
            c = c + 1;
        end
        boxes{l,1} = imbinarize(b);
    end
%{
hold on;
for k=1:size(boxes,1)
    curr = boxes{k,1};
    [L,~] = bwlabel(curr);
    box = regionprops(L, 'BoundingBox');
    rectangle('Position',box(1).BoundingBox,'EdgeColor','yellow');
end
hold off;

    
    temp = imshow(image); axis on;
    hold on;
    [L,~] = bwlabel(boxes{1});
    box = regionprops(L, 'BoundingBox');
    rectangle('Position',box(1).BoundingBox,'EdgeColor','yellow');
    hold off;
    saveas(temp,strcat('..\',sprintf('%d',i+1)),'png');
    %}
    
    for j=1:size(windowcms,1)
        lab = rgb2lab(image(windowindeces{j,1}(1,1):windowindeces{j,2}(1,1), windowindeces{j,1}(1,2):windowindeces{j,2}(1,2),:));
        l = lab(:,:,1);
        a = lab(:,:,2);
        b = lab(:,:,3);
        
        gmmfg = fitgmdist([double(l(shapemodels{j,1} > .75)) double(a(shapemodels{j,1} > .75)) double(b(shapemodels{j,1} > .75))],3,'Options',statset('MaxIter',7500),'RegularizationValue',.00000000001);
        gmmbg = fitgmdist([double(l(shapemodels{j,1} > .25)) double(a(shapemodels{j,1} > .25)) double(b(shapemodels{j,1} > .25))],3,'Options',statset('MaxIter',7500),'RegularizationValue',.00000000001);
        
        l = double(reshape(l,[size(l,1)*size(l,2),1]));
        a = double(reshape(a,[size(l,1)*size(l,2),1]));
        b = double(reshape(b,[size(l,1)*size(l,2),1]));
        
        pdffg = pdf(gmmfg,[l a b]);
        pdfbg = pdf(gmmbg,[l a b]);
        
        avgF = (pdffg + pdf(oldGFB{j,1}, [l a b]))./2;
        avgB = (pdfbg + pdf(oldGFB{j,2}, [l a b]))./2;
        
        newcm = avgF./(avgF + avgB);
        newcm = reshape(newcm,[size(windows{j,1},1) size(windows{j,1},2)]);
   
        if size(newcm(newcm > 0),1) < size(windowcms{j,1}(windowcms{j,1} > 0),1)
            windowcms{j,1} = newcm;
            localfg = foreground(windowindeces{j,1}(1,1):windowindeces{j,2}(1,1), windowindeces{j,1}(1,2):windowindeces{j,2}(1,2));
            localdist = disttrans(windowindeces{j,1}(1,1):windowindeces{j,2}(1,1), windowindeces{j,1}(1,2):windowindeces{j,2}(1,2));
            numerator = 0;
            denominator = 0;
            for r=1:size(windowcms{j,1},1)
                for c=1:size(windowcms{j,1},2)
                    L = localfg(r,c);
                    p = windowcms{j,1}(r,c);
                    d = localdist(r,c);
                    w = exp(-(d^2)/(round(size(windowcms{j,1},1)/2)));
                    numerator = numerator + abs(L - p) * w;
                    denominator = denominator + w;
                end
            end
            windowcolorconf{j,1} = 1 - (numerator/denominator);
            oldGFB{j,1} = gmmfg;
            oldGFB{j,2} = gmmbg;
        end
    end
    %{
    temp = imshow(windowcms{1});
    saveas(temp,strcat('..\',sprintf('%d',7+i)),'png');
    %}
    
    view = imref2d([size(images{i},1) size(images{i},2)]);
    foreground = imwarp(foreground,tform,'OutputView',view);
    outline = bwmorph(foreground,'remove');
    disttrans = bwdist(outline);
    
    %{
    temp = imshow(foreground(windowindeces{1,1}(1,1):windowindeces{1,2}(1,1), windowindeces{1,1}(1,2):windowindeces{1,2}(1,2)));
    saveas(temp,strcat('..\',sprintf('%d',13+i)),'png');
    %}
    
    for j=1:size(windowindeces,1)
        shapemodels{j,1} = zeros(size(windows{j,1},1), size(windows{j,1},2));
        for r=windowindeces{j,1}(1,1):windowindeces{j,2}(1,1)
            for c=windowindeces{j,1}(1,2):windowindeces{j,2}(1,2)
                L = foreground(r,c);
                d = -(disttrans(r,c))^2;
                a = (size(windowcms{j,1},1) - 2)/((1 - .85)^2);
                if windowcolorconf{j,1} > .85
                    s = 2 + (a * (windowcolorconf{j,1} - .85)^2);
                else
                    s = 2;
                end
                shapemodels{j,1}(r-windowindeces{j,1}(1,1)+1,c-windowindeces{j,1}(1,2)+1) = 1 - exp(d/(s^2));
            end
        end
    end
    %{
    temp = imshow(shapemodels{1});
    saveas(temp,strcat('..\',sprintf('%d',19+i)),'png');
    %}
    
    fgmaps = cell(size(windowindeces,1),1);
    for j=1:numwindows
        L = foreground(windowindeces{j,1}(1,1):windowindeces{j,2}(1,1), windowindeces{j,1}(1,2):windowindeces{j,2}(1,2));
        fgmaps{j,1} = (shapemodels{j,1}.*L) + ((1 - shapemodels{j,1}).* windowcms{j,1});
    end
%{
    temp = imshow(fgmaps{1});
    saveas(temp,strcat('..\',sprintf('%d',24+i)),'png');
%}
    
    for j=1:numwindows
       whole = zeros(size(images{i},1), size(images{i},2));
       for r=windowindeces{j,1}(1,1):windowindeces{j,2}(1,1)
           for c=windowindeces{j,1}(1,2):windowindeces{j,2}(1,2)
               whole(r,c) = fgmaps{j,1}(r-windowindeces{j,1}(1,1)+1,c-windowindeces{j,1}(1,2)+1); 
           end
       end
       fgmaps{j,1} = whole;
    end
    
    finalfgmask = zeros(size(images{i},1), size(images{i},2));
    for r=1:size(finalfgmask,1)
        for c=1:size(finalfgmask,2)
            numerator = 0;
            denominator = 0;
            for k=1:numwindows   
                if r >= windowindeces{k,1}(1,1) && r <= windowindeces{k,2}(1,1) && c >= windowindeces{k,1}(1,2) && c <= windowindeces{k,2}(1,2)
                    center = windowindeces{k,3};
                    distance = sqrt((c-center(1,2))^2 + (r-center(1,1))^2);
                    numerator = numerator + (fgmaps{k,1}(r,c) * (distance + .1)^(-1));
                    denominator = denominator + ((distance + .1)^(-1));
                end
            end
            finalfgmask(r,c) = numerator/denominator;
        end
    end
    %{
    temp = uint8(finalfgmask.*255);
    temp = imshow(temp);
    saveas(temp,strcat('..\Output\',sprintf('%d',i)),'png');
    %}
      
    finalfgmask = imfill(finalfgmask > 0,'holes');
    lazybw = lazysnapping(image,superpixels(image,1000),finalfgmask,~finalfgmask);
    finalfgmask = bwmorph(lazybw,'remove');
    
    temp = image;
    R = temp(:,:,1);
    G = temp(:,:,2);
    B = temp(:,:,3);
    R(finalfgmask) = 255;
    G(finalfgmask) = 0;
    B(finalfgmask) = 0;
    temp(:,:,1) = R;
    temp(:,:,2) = G;
    temp(:,:,3) = B;
    temp = imshow(temp);
    saveas(temp,strcat('..\Output\',sprintf('%d',i)),'png');
end




