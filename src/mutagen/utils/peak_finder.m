
function partials_list = peak_finder(spectrum, mwidth)
    
    spectrum = spectrum(:);
    
    % Enforce odd mwidth
    mwidth = floor(mwidth/2) * 2 + 1;

    mag_spectrum = abs(fft(padded_in));
    half_spectrum = mag_spectrum(1:end/2);
    half_spectrum_size = size(half_spectrum, 1);
    
    moving_window = zeros(half_spectrum_size-mwidth+1, mwidth);
    for k = 1:mwidth
        moving_window(:,k) = half_spectrum(k:end-mwidth+k, 1);
    endfor
    
    med_filt_output = median(moving_window, 2);
    % Many values in salient_values will be zero due to the subtraction here.
    salient_values_start = (mwidth+1)/2;
    salient_values_end = -salient_values_start + 1;
    salient_values = ...
        log(half_spectrum(salient_values_start:end+salient_values_end)+1e-300) - ...
        log(med_filt_output+1e-300);
    % Because of all the zeros in salient_values, we need to deliberately reject peaks
    % at height zero, hence the third factor in this multiplication.
    indicators_of_pos_peaks = (salient_values(3:end) <= salient_values(2:end-1)) .* ...
                              (salient_values(2:end-1) > salient_values(1:end-2)) .* ...
                              (salient_values(2:end-1) > 0);
    peak_indices = find(indicators_of_pos_peaks);
    [sorted_peak_heights peak_perm] = sort(salient_values(peak_indices));
    accum_energy_by_freq = cumsum(sorted_peak_heights .* sorted_peak_heights);
    normalized_accum = sqrt(accum_energy_by_freq/accum_energy_by_freq(end)) * ...
                       half_spectrum_size;
    % We're trying to find a few peaks that account for more than their share of the energy
    % in the spectrum. normalized_accum is a non-decreasing series with a non-decreasing
    % derivative that captures how flat the distribution is. (It's something like the
    % integral of a histogram.) A purely flat spectrum will result in a flat diagonal. The
    % more "peaky" the distribution, the more prominent an elbow there will be. Our ad hoc
    % criteria are to say that the elbow is the point where the slope exceeds 2, and then to
    % take points beyond the elbow that satisfy y <= 4x-3 when the axis intervals are
    % normalized to a square with corners {0|1},{0|1}.
    % It probably is possible to do this with the derivative of this function, but we're not
    % going to put the time in for that right now, as this approach already works.
    crook_of_elbow = find((normalized_accum(3:end) - normalized_accum(1:end-2)) > 4)(1) + 1;
    % y/half_spectrum_size <= 4(x/half_spectrum_size) - 3
    % y <= 4x - 3*half_spectrum_size
    first_point = find(normalized_accum(crook_of_elbow:end) <= ...
                       4*[crook_of_elbow:half_spectrum_size]' - 3*half_spectrum_size)(1) + ...
                  crook_of_elbow - 1;   
keyboard();

    % The sorting permutation tells us which bins in the original spectrum these peaks come
    % from.
    strongest_peak_bins = peak_perm(first_point:end);
    % Adjust these indices, which point into salient_values, so they point properly into
    % half_mag_spectrum.
    strongest_peak_bins = strongest_peak_bins + salient_values_start - 1;
    c0 = half_spectrum(strongest_peak_bins);
    to_the_left = half_spectrum(strongest_peak_bins-1);
    to_the_right = half_spectrum(strongest_peak_bins+1);
    c1 = 0.5 * (to_the_right - to_the_left);
    c2 = 0.5 * (to_the_right + to_the_left) - c0;
    xs_at_peaks = -c1 ./ (2 * c2);
    peak_strengths = (c2 .* xs_at_peaks + c1) .* xs_at_peaks + c0;
    % -1: Indices need to be zero-based, not one-based, because the first bin corresponds
    %     to 0 Hz.
    %partials_list = [4096/81*8 1024; 4096/64*8 1024; 4096/49*8 1024];
    partials_list = [xs_at_peaks+strongest_peak_bins-1 peak_strengths];
    
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = peaks_and_troughs(in)

    in = in(:);
    
    first_diff = diff(in);
    signs_of_diff = sign(first_diff);
    % A diff value with zero sign takes on the subsequent sign.
    zero_signs = find(signs_of_diff(1:end-1) == 0);
    signs_of_diff(zero_signs) = signs_of_diff(zero_signs + 1);
    % Positive values of diff_zero_crossings indicate peaks; negatives are troughs.
    diff_zero_crossings = find(sign(diff(signs_of_diff)));
    
    % UNFINISHED

endfunction
