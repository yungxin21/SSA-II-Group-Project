function corrected_img = limbDarkness(fits_input,plot_spacing,corrected_fig_num)
    img = fitsread(fits_input, "image");

    % Magnetogram option (in separate branch)

    % AIA EUV option
    % euv_lev1 = fitsread(fits_input, "image");
    % figure(3); imagesc(euv_lev1); colormap("bone"); colorbar;

    % 1. Identify the solar discs

    info = fitsinfo(fits_input);
    img_keywords = info.Image.Keywords;
    valIdx = strcmp(img_keywords(:,1), 'CRPIX1');
    valIdx2 = strcmp(img_keywords(:,1), 'CRPIX2');
    center_x = round(img_keywords{valIdx, 2});
    center_y = round(img_keywords{valIdx2, 2});
    
    % find the first and last index with positive values larger than 1000
    xmin=find(img(center_x,:)/1000>=1,1,"first"); 
    xmax=find(img(center_x,:)/1000>=1,1,"last");
    diameter_x=xmax-xmin+1;
    radius_x=diameter_x/2;
    
    ymin=find(img(:,center_y)/1000>=1,1,"first");
    ymax=find(img(:,center_y)/1000>=1,1,"last");
    diameter_y=ymax-ymin+1;
    radius_y=diameter_y/2;
    
    radius = floor(mean([radius_x,radius_y])); % the radius of the sun in pixels
    
    % 2. Correct for limb darkening.
    % SDO HMI continuum image wavelength: 617.3 nm (limb darkness is
    % dependent on captured wavelength)
    
    % grid of pixel coordinates
    [X, Y] = meshgrid(1:size(img, 2), 1:size(img, 1));
    % calculate the radial distance (r) from the center for every pixel
    r = sqrt((X - center_x).^2 + (Y - center_y).^2);
    
    % here mu = cos(theta), where theta is the heliocentric angle
    mu = sqrt(max(0, 1 - (r ./ radius).^2)); 
    % limb darkening coefficients (see references):
    u1 = 0.836; 
    u2 = -0.204;
    
    % equals ~1 at the around the centre (i.e, no correction) and near 0 around
    % the edges (a lot of correction)
    limb_fit = 1 - u1 - u2 + u1 * mu + u2 * (mu).^2;
    
    % need to change this to work with sunspot.m
    solar_disk_mask = (r <= radius); % basically just 0 where the sun isnt
    % and 1 where the sun is, so correction is only applied on sun and not background
    % solar_disk_mask(solar_disk_mask<0) = NaN;
    
    % apply the corrected intensity values to the array of 0s and 1s
    corrected_img = img;
    corrected_img(solar_disk_mask) = img(solar_disk_mask) ./ limb_fit(solar_disk_mask);

    % block of code to plot corrected image
    if (corrected_fig_num > 0) 
        figure(corrected_fig_num);
        imagesc(corrected_img)
        set(gca, "YDir", "normal");
        colormap("bone");
        colorbar;
    end

    % adjustable intensity plot with normalized radii.
    if (plot_spacing > 0)
        correctedModes = [];
        originalModes = [];
        % view intensity over r/r_solar plot.
        for i=1:plot_spacing:radius
            ring = (round(r)==round(i));
            correctedModes = [correctedModes, mode(corrected_img(ring))]; 
            originalModes = [originalModes, mode(img(ring))];
        end
        figure(2); plot((1:plot_spacing:radius)/radius,correctedModes); hold on;
        plot((1:plot_spacing:radius)/radius,originalModes); hold off;
        title('Intensity over Radius Plot')
        xlabel('r/r_☉')
        ylabel('intensity')
        xlim([0 1.25])
        ylim([0 55000])
    end

    % retaining HDU in FITS files can only be achieved using low-level fits
    % code
    % checks if fits_output exists, then deletes it.
    % if exist(fits_output, 'file') == 2
    %     delete(fits_output);
    % end
    
    % % Initialize the Primary Image HDU matching data type and size
    % % fptr = fits.openFile("corrected.fits","readwrite");
    % fptr = fits.createFile(convertStringsToChars(fits_output)); 
    % 
    % imgSize = size(corrected_img);
    % fits.createImg(fptr, 'double_img', imgSize);
    % fits.writeImg(fptr, corrected_img);
    % 
    % % Filter out automatic structural headers (SIMPLE, BITPIX, NAXIS) to prevent duplicates
    % trueKeys = {'WCSAXES', 'MJD-OBS', 'DATE-OBS', 'CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'};
    % for i = 1:size(img_keywords, 1)
    %     img_keyName = img_keywords{i, 1};
    %     img_keyValue = img_keywords{i, 2};
    %     img_keyComment = img_keywords{i, 3};
    % 
    %     % Only write custom metadata tags, avoiding reserved structural keys
    %     if ~isempty(img_keyName) && ismember(img_keyName, trueKeys)
    %         if islogical(img_keyValue)
    %             if img_keyValue
    %                 img_keyValue = 'T';
    %             else 
    %                 img_keyValue = 'F';
    %             end
    %         end
    %         if isempty(img_keyValue), img_keyValue = ' '; end
    %         if isempty(img_keyComment), img_keyComment = ' '; end
    %         if isstring(img_keyValue), img_keyValue = char(img_keyValue); end
    %         if isstring(img_keyComment), img_keyComment = char(img_keyComment); end    
    %         try
    %             fits.writeKey(fptr, img_keyName, img_keyValue, img_keyComment);
    %         catch ME
    %             warning('Failed to write key "%s" with value class "%s": %s', ...
    %                 img_keyName, class(img_keyValue), ME.message);
    %         end
    %     end
    % end
    % fits.closeFile(fptr)
end