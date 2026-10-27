% changes to be made ON LOCAL DEVICE
% 1. In each unzipped folder, there is a "index.txt" file. 
% Go in it, and delete the lines from the start until "#DATA".
% 2. In the lines for intensitylist and magnetoList, change the folder to
% your dataset folder, then "\index.txt", like how ive done.
% 3. uncomment the for loop lines
% 4. in the limbCorrection function line (first line after where the for loop starts), 
% change the long "...fits" file to intensityList(i).
% Also change the input in the "magimg=..." line to magnetoList(i)

% read .txt file lines and return an array with all the file names
% intensityList = readlines("JSOC_20261024_000868\index.txt");
% magnetoList = readlines("JSOC_20261025_004772\index.txt");
% initialization to prevent error for the first concatenation
spotInfoList = [];

% temp variables
intensityTemp = "hmi.Ic_720s.20240504_040000_TAI.3.continuum.fits";
magTemp = "hmi.M_720s.20240504_040000_TAI.3.magnetogram.fits";

% for loop for automation
% for i=1:size(file_names,1)
info = fitsinfo(intensityTemp); % change this line
keywords = info.Image.Keywords; % retrieve info from HDU file.
valIdx1 = strcmp(keywords(:,1), 'CDELT1');
valIdx2 = strcmp(keywords(:,1), 'CDELT2');
pixelScale1 = keywords{valIdx1, 2}; 
pixelScale2 = keywords{valIdx2, 2};
avgPxScale = (pixelScale1+pixelScale2)/2;

[img, r_sun, u_center, v_center] = limbCorrection(intensityTemp,0,0);
magimg = fitsread(magTemp, "image");
tfimg = tf2square(img, r_sun, u_center, v_center, 4096, 4096);
% figure(6); imagesc(tfimg); colormap; colorbar; set(gca, 'YDir', 'normal');
[individuals, GE] = sunspotClass(tfimg, r_sun, u_center, v_center, avgPxScale); % intensitygram

magimg = tf2square(magimg, r_sun, u_center, v_center, 4096, 4096);

% [individuals, GE, clippedImg, penumbraMask, umbraMask] = sunspotClass(img, r_sun, u_center, v_center, avgPxScale); % intensitygram
[magIndivs, magGE] = haleslaw(tfimg, magimg, r_sun, u_center, v_center, avgPxScale); % magnetogram

% Wolf Number
% will be used in the final main function. f and g are outputs (or just
% Wolf Number).

f1 = numel(individuals);
g1 = numel(GE);
f2 = numel(magIndivs);
g2 = numel(magGE);

sunspotNum1 = 10*g1 + f1;
sunspotNum2 = 10*g2 + f2;

% not storing image info for now, maybe future implementation.
intSpotInfo = struct('sunspotNum', sunspotNum1, 'individualData', individuals, 'groupData', GE); % create struct to store sunspot information
magSpotInfo = struct('sunspotNum', sunspotNum2, 'individualData', magIndivs, 'groupData', magGE);
spotInfoList = [spotInfoList; intSpotInfo, magSpotInfo]; % concatenation

% could also create a class with no methods and input that into a
% struct to define fields, however it requires defining datatypes of each field

% end