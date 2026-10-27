function [maskedImg, mask] = sigmaClip(data, shadowID, maxI)
    % Make all negative values to be NaN for ease of finding 0 in data
    maskedImg = double(data);
    maskedImg(maskedImg<=0) = NaN;
    % Dedicated map of where NaN's are
    nanMask = isnan(maskedImg);
    % Start mask (true = valid, false = clipped)
    mask = true(size(maskedImg));
    for i=1:maxI
        % Make current unclipped values as 1D vector
        validData = maskedImg(mask & ~nanMask);
        if isempty(validData), break; end
        % Running statistics on valid data and omitting clipped values
        bins = 0:0.00001:1.3;
        [count,bins] = hist(validData, bins);
        [~, maxIdx] = max(count);
        % disp("maxVal is " + maxVal);
        % disp("max index is " + maxIdx);
        centreVal = bins(maxIdx);
        % stdVal = std(validData);

        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % validation shows error. Needs recoding
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % lower bound? of possible quantiles
        % pen = 0.001, 0.999
        % umb = 0.00025, 0.99975
        if shadowID == 1 % this is the penumbra case (i.e. shadowID == 1)
            quantileVal = 0.015; % still hardcoded
        elseif shadowID == 2 % this is the umbra case (shadowID == 2)
            quantileVal = 0.001;
        else % any other shadowID, if entered in error or on purpose.
            lowBound = 0;
            upBound = 1.3;
        end
        if shadowID > 0
            lowBound = quantile(validData, quantileVal);
            upBound = quantile(validData, 1-quantileVal);
        end
        FW = upBound - lowBound;

        % % debugging
        % fprintf('n%% Width (Empirical Quantile): %.4f\n', FW);
        % fprintf('Interval centered at %g: [%.4f, %.4f]\n', median(validData), lower_bound, upper_bound);

        % % back when standard deviation was used
        % disp(stdVal);
        % % find nSigma using HWHM?
        % fracMax = maxVal / 100;
        % peakIndices = find(count >= fracMax);
        % nSigma = bins(peakIndices(end))-bins(peakIndices(1));
        % disp(nSigma);

        % if stdVal == 0, break; end
        % % disp("isrunning2");
        % Defining boundaries
        lBound = centreVal - FW;
        % Use bounds to generate new mask 
        newMask = mask & (maskedImg >= lBound);
        if isequal(newMask,mask), break; end
        % disp("isrunning3");
        mask = newMask;
    end
    % Replacing clipped values to keep same matrix dimensions
    maskedImg(~mask & ~nanMask) = -0.0923;
end