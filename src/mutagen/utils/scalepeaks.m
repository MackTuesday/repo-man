
function out = scalepeaks(spectrum, mwidth, threshold, scale)
    
    spectrum = spectrum(:);
    magspectrum = abs(spectrum);

    mmidpoint = (mwidth-1) / 2 + 1;
    ex_half_spectrum = [magspectrum(mmidpoint:-1:2); magspectrum(1:end/2+mmidpoint)];
    ex_half_spectrum_size = size(ex_half_spectrum, 1);
    
    moving_window = zeros(ex_half_spectrum_size-mwidth+1, mwidth);
    for k = 1:mwidth
        moving_window(:,k) = ex_half_spectrum(k:end-mwidth+k, 1);
    endfor
    
    med_filt_output = median(moving_window, 2);
    
    % Many values in salient_values will be zero due to the subtraction here.
    salient_values_start = (mwidth+1)/2;
    salient_values_end = -salient_values_start + 1;
    salient_values = ...
        log(ex_half_spectrum(salient_values_start:end+salient_values_end)+1e-300) - ...
        log(med_filt_output+1e-300);
        
    % Because of all the zeros in salient_values, we need to deliberately reject peaks
    % at height zero, hence the third factor in this multiplication.
    % indicators_of_pos_peaks will be a bunch of boolean values, and peak_indices will point
    % into salient_values
    indicators_of_pos_peaks = (salient_values(3:end) <= salient_values(2:end-1)) .* ...
                              (salient_values(2:end-1) > salient_values(1:end-2)) .* ...
                              (salient_values(2:end-1) > 0);
    peak_indices = find(indicators_of_pos_peaks) + 1;

    % Apply quadratic interpolation
    c0 = salient_values(peak_indices);
    to_the_left = salient_values(peak_indices-1);
    to_the_right = salient_values(peak_indices+1);
    c1 = 0.5 * (to_the_right - to_the_left);
    c2 = 0.5 * (to_the_right + to_the_left) - c0;
    xs_at_peaks = -c1 ./ (2 * c2);
    peak_strengths = (c2 .* xs_at_peaks + c1) .* xs_at_peaks + c0;
    
    % Identify peaks that are strong enough; these indices point into peak_strengths and
    % therefore also into xs_at_peaks and peak_indices.
    indices_of_strong_enough_peaks = find(peak_strengths > threshold);
     
    % For each peak, the two adjacent bins that contribute to that peak; both xs_at_peaks
    % and peak_indices point into salient_values.
    exact_peak_indices = xs_at_peaks(indices_of_strong_enough_peaks) + ...
                         peak_indices(indices_of_strong_enough_peaks);
                         
    % Some values of adjacent_bins_left will be equal to corresponding values of
    % adjacent_bins_right.
    adjacent_bins_left = ...
        min(floor(exact_peak_indices), peak_indices(indices_of_strong_enough_peaks));
    adjacent_bins_right = ...
        max(peak_indices(indices_of_strong_enough_peaks), ceil(exact_peak_indices));
    
    % More bins surrounding the peak
    bins_to_scale = union(adjacent_bins_left, adjacent_bins_right);
    for i = 1:mmidpoint-1
        bins_to_scale = union(bins_to_scale, adjacent_bins_left-i);
        bins_to_scale = union(bins_to_scale, adjacent_bins_right+i);
    endfor

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

endfunction
