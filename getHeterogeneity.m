function [dataMat, numDoublePos, numRedOnly, numGreenOnly, THRESHGREEN]...
    = getHeterogeneity(bf, red, green)

    close all;
    
    Ibf = imread(bf);
    %%%imtool(Iphase,[]);
    Ired = imread(red);
    %%%imtool(Ired,[]);
    Igreen = imread(green);
    %%%imtool(Igreen,[]);
    
%     greenref = imread('AN1 CF fixed GREEN 1s 9 this one.tif');
%     greenrefd = im2double(greenref);
%     greenrefg = rgb2gray(greenrefd);
    
    %merge = imfuse(Ired, Igreen, 'blend', 'scaling', 'joint');
    %%%imtool(merge,[]); %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    %%convert microscopy images to doubles
    Ibfd = im2double(Ibf);
    Iredd = im2double(Ired);
    Igreend = im2double(Igreen);
    
    %%convert double images to grayscale
    Ibfg = imcomplement(mat2gray(Ibfd));
    Iredg = mat2gray(Iredd);
    Igreeng = mat2gray(Igreend);

    %%calcuate threshold values for each color
    THRESHBF = graythresh(Ibfg);    
    THRESHRED = graythresh(Iredg);
    THRESHGREEN = graythresh(Igreeng); %graythresh(Igreeng);

    %%set stringency of calculated threshold (>1: more stringent, <1: less)
    BFSTRICT = 1;
    REDSTRICT = 1.2;
    GREENSTRICT = 2;

    %%generates masked bw images
    if THRESHBF*BFSTRICT < 1
        bwbf = im2bw(Ibfg, THRESHBF*BFSTRICT);
    else bwbf = im2bw(Ibfg, THRESHBF);
    end
    if THRESHRED*REDSTRICT < 1
        bwred = im2bw(Iredg, THRESHRED*REDSTRICT);
    else bwred = im2bw(Iredg, THRESHRED);
    end
    if THRESHGREEN*GREENSTRICT < 1
        bwgreen = im2bw(Igreeng, THRESHGREEN*GREENSTRICT);
    else bwgreen = im2bw(Igreeng, THRESHGREEN);
    end
    
    %imtool(bwbf,[]);
    %imtool(bwred,[]);
    %imtool(bwgreen,[]);

    %%labels each region with specific pixel values
    bflabel = bwlabel(bwbf);
    redlabel = bwlabel(bwred);
    greenlabel = bwlabel(bwgreen);

    for i = 1:1040
        for j = 1:1392
            if bflabel(i, j) == 0
                Ibfg(i, j) = 0;
                Iredg(i, j) = 0;
                Igreeng(i, j) = 0;
            end
        end
    end
    
%     %imtool(Ibfg,[]);
%     %imtool(Iredg,[]);
%     %imtool(Igreeng,[]);

    cellprops = regionprops(bwbf, 'Area', 'PixelIdxList');
    cellpropsize = size(cellprops);
    numcellclumps = cellpropsize(1);
    singlecells = bflabel;
    %%imtool(singlecells,[]);
    
    imsize = size(singlecells);
    x = imsize(1);
    y = imsize(2);

    %%eliminates cell clumps and relabels
    for i = 1:numcellclumps
        if cellprops(i).Area > 200 %empirical value **********************
            singlecells(singlecells==i)=0;
        elseif cellprops(i).Area < 20 %same
            singlecells(singlecells==i)=0;
        end
    end
    
    %%gets rid of red and green light that belongs to clumps
    for i = 1:x
        for j = 1:y
            if singlecells(i, j)==0
                greenlabel(i, j)=0;
                redlabel(i, j)=0;
            end
        end
    end

    singlecells = bwlabel(imclearborder(singlecells));
    singleprops = regionprops(singlecells, 'PixelIdxList', 'Area');
    cellnumsize = size(singleprops);
    numcells = cellnumsize(1);
    %%imtool(singlecells,[]);
    
    %%imtool(greenlabel,[]);

    greenspots{numcells, 1} = [];
    redspots{numcells, 1} = [];

    %%finds red/green regions within each cell and binds to cellnum
    for i = 1:numcells
        location = singleprops(i).PixelIdxList;
        uniquegreen = unique(greenlabel(location));
        uniquegreen(uniquegreen == 0) = [];
        greenspots(i, 1) = {uniquegreen'};
        uniquered = unique(redlabel(location));
        uniquered(uniquered == 0) = [];
        redspots(i, 1) = {uniquered'};
    end
    
    greenspots = greenspots(~cellfun('isempty',greenspots));
    redspots = redspots(~cellfun('isempty',redspots));
    
    dimensiongreen = size(greenspots);
    dimensionred = size(redspots);
    sizegreen = dimensiongreen(1);
    sizered = dimensionred(1);
    
    numDoublePos = sizegreen;
    numGreenOnly = 0;

    %%initializes a matrix for use in PCA each row is one cell's data
    dataMat = zeros(numcells, 5); %change cols for number of params
    %area 1
    %max GFP intensity 2
    %mean GFP intensity 3
    %max RFP intensity 4
    %mean RFP intensity 5
    % Obviously more parameters could be mined here, just add them to the
    % following for loop and change the initial matrix above

    %%each ROW is parameter values for a single cell
    %%each COLUMN is a single parameter's values for all cells
    %%matrix is for ONE image, append matrices vert. from mult conditions
    for i = 1:numcells
        dataMat(i, 1) = singleprops(i).Area; % cell area
        location = singleprops(i).PixelIdxList;
        greengrab = Igreend(location);
        redgrab = Iredd(location);
        dataMat(i, 2) = max(greengrab);
        dataMat(i, 3) = mean(greengrab);
        dataMat(i, 4) = max(redgrab);
        dataMat(i, 5) = mean(redgrab);
    end
    
    for i = 1:numcells
        location = singleprops(i).PixelIdxList;
        greenbinary = bwgreen(location);
        redbinary = bwred(location);
        if mean(redbinary)== 0
            if mean(greenbinary) ~= 0
                numDoublePos = numDoublePos - 1;
                numGreenOnly = numGreenOnly + 1;
            end
        end
    end
    
    numRedOnly = sizered-numDoublePos;
    
end