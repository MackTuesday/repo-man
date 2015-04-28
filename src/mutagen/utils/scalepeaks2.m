# This file has been abandoned. It contains no usable code.

# function troughindices = trough_analysis(X)
function troughindices = scalepeaks2(X, order, threshold)

    X = X(:);
    sizex = size(X, 1);

end


#    levelnumber = 1;
#    troughindices{levelnumber} = [1:sizex]';
#   levelsize = sizex;
#    while (levelsize > 4 || (levelsize == 4 && (level(2) != level(3))) || ...
#           (levelsize == 3 && (level(2) >= level(1) || level(2) > level(3))))
#        level = X(troughindices{levelnumber});
#        levelnumber = levelnumber + 1;
#        levelsize = size(level, 1)
#        tempindices = [1; find((level(1:end-2) > level(2:end-1) & ...
#                                level(2:end-1) <= level(3:end)) | ...
#                               (level(1:end-2) >= level(2:end-1) & ...
#                                level(2:end-1) < level(3:end))) + 1; levelsize];
#        troughindices{levelnumber} = troughindices{levelnumber-1}(tempindices);
#    end

#    for i = 1:levelnumber
#        level = X(troughindices{i});
#        level(1) = min(level(1), level(2));
#        level(end) = min(level(end), level(end-1));
#    end

#    areas = peak_areas(X, troughindices, sizex/4);
    
#keyboard();
#end

function areas = peak_areas(X, troughindices, resolution)

    X = X(:);
    sizex = size(X, 1);

    troughindices = troughindices(:);
    troughlevels = size(troughindices, 1);
    for i = 1:troughlevels
        intervals = troughindices{i}(2:end) - troughindices{i}(1:end-1);
        if (max(intervals) < resolution)
            numlevels = i;
        else
            break;
        end
    end
    
    X = X - min(X);
interpolation = zeros(sizex, 1);
    for i = numlevels-1:-1:1
        level = X(troughindices{i});
        interpolation = interp1(troughindices{i}, level, [1:sizex]');
        trapareas = 0.5 * (interpolation(1:end-1) + interpolation(2:end));
        cumarea = [0 ; cumsum(trapareas)];
        areas{i} = cumarea(troughindices{i}(2:end)) - cumarea(troughindices{i}(1:end-1));
figure(1);
plot(X);
        X = X - interpolation;
figure(2);
plot(X);       
keyboard();
    end    
    
end

#    interps = zeros(sizex, levelnumber);
#    for i = 1:levelnumber
#        level = X(troughindices{i});
#        level(1) = min(level(1), level(2));
#        level(end) = min(level(end), level(end-1));
#        interps(:,i) = interp1(troughindices{i}, level, [1:sizex]');
#    end


function out = scalepeaksx(spectrum, mwidth, threshold, scale)
    
    spectrum = spectrum(:);
    magspectrum = abs(spectrum);

    mmidpoint = (mwidth-1) / 2 + 1;
    ex_half_spectrum = [magspectrum(mmidpoint:-1:2); magspectrum(1:end/2+mmidpoint)];
    ex_half_spectrum_size = size(ex_half_spectrum, 1);
    
    [hierarchy sizes] = trough_hierarchy(ex_half_spectrum);
    
    numlevels = size(sizes, 1);
    if (numlevels < 2)
        out = spectrum;
        return;
    end

    % For each level 1:end-1
    %   Get sums of bins between troughs (heights are from interpolated inter-trough height, 
    %    not from zero)
    %   Get average of all sums between troughs at the next level of the hierarchy
    %   Suppress bins strictly between troughs where corresponding sum exceeds the threshold
    troughindices = hierarchy(1:sizes(1));
    troughspacing = diff(troughindices);
    for level = 2:numlevels
    endfor
    
    % Some values of adjacent_bins_left will be equal to corresponding values of
    % adjacent_bins_right.
    adjacent_bins_left = ...
        min(floor(exact_peak_indices), peak_indices(indices_of_strong_enough_peaks));
    adjacent_bins_right = ...
        max(peak_indices(indices_of_strong_enough_peaks), ceil(exact_peak_indices));
    
    % More bins surrounding the peak
    bins_to_scale = [adjacent_bins_left; adjacent_bins_right];
    for i = 1:mmidpoint-1
        bins_to_scale = [bins_to_scale; adjacent_bins_left-i];
        bins_to_scale = [bins_to_scale; adjacent_bins_right+i];
    end

    bins_to_scale = min(bins_to_scale, size(salient_values, 1));
    bins_to_scale = max(1, bins_to_scale);
    bins_to_scale = unique(bins_to_scale);
    
    % These "adjacent_bins_*" arrays are indices that point into salient_values, therefore
    % they also point into med_filt_output. This means they also point properly into
    % spectrum.
    spectrum(bins_to_scale) = ...
        med_filt_output(bins_to_scale) .* exp(scale * salient_values(bins_to_scale)) .* ...
        spectrum(bins_to_scale) ./ abs(spectrum(bins_to_scale));
    spectrum(end/2+2:end) = conj(spectrum(end/2:-1:2));
    
    out = spectrum;

end

