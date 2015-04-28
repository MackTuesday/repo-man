
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grains, residual] = MPvR(in, window_length)

    in = in(:);
    input_size = size(in,1);
    
    % Pad input with half a window length at beginning maybe more at the end
    padded_input_size = ceil((input_size + window_length) / window_length) * window_length;
    %padded_input = [ zeros(window_length/2,1); in; ...
    %                 zeros(padded_input_size - input_size - window_length/2, 1) ];
    padded_input = [ in; ...
                     zeros(padded_input_size - input_size - window_length/2, 1) ];
    
    % For each half window length
    stride_size = window_length;
    num_strides = padded_input_size / stride_size;
    for stride = 0:num_strides-1
        % Section goes from here to here + window length
        section_begin = stride * stride_size + 1;
        section_end = section_begin + stride_size - 1;
        section = padded_input(section_begin:section_end);
        
        residual = section;
        
        oversample_exp = 3;
        
        % Do the following until there don't seem to be any salient partials
        %while (0 == 0)
        
            % Take windowed FFT
            window_function = 0.5 * (1 + cos([-window_length/2+1:window_length/2]' * ...
                                             2 * pi / window_length));
            windowed_section = residual .* window_function;
            spectrum = superfft(windowed_section, oversample_exp);

            % Find salient partials
            target_count = 3;
            partials_list = find_partials(spectrum, target_count);
            
            % Build basis vectors from all partials
            max_grain_length = window_length / 10;
            time_resolution = 100;
            [basis whole_grain_indices] = ...
                build_basis(residual, spectrum, partials_list, window_length, ...
                            max_grain_length, time_resolution);
            
            ls_basis = basis + 1e9;
            ls_basis = ls_basis - 1e9;
            
            % Do least squares
            coeffs = ridge_regression(ls_basis, residual, 0.001);

            % Build approximation
            appx = basis(:,whole_grain_indices) * coeffs(whole_grain_indices,1);
keyboard();
            
            % Subtract to get residual
            residual = residual - appx;
            
        %endwhile
        
    endfor

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = superfft(in, oversample_exp)
    
    in = in(:);

    if (oversample_exp < 1)
        out = fft(in);
        return;
    endif
    
    in_length = size(in, 1);
    oversample_factor = 2^oversample_exp;
    padded_in = [in ; zeros((oversample_factor-1)*in_length, 1)];
    out = fft(padded_in);
    
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function partials_list = find_partials(spectrum, target_count)

    partials_list = [];
    
    threshold = 0.3;
    peaks_found = 0;
    threshold_step = 1 / sqrt(2);
    mag_spectrum = abs(spectrum);

    while (threshold > 0.01 && peaks_found < target_count)
        [peaks holes] = RED(log(mag_spectrum(2:floor(end/2)-1)'+1e-300), threshold);
        peaks = peaks';
        peaks_found = size(peaks, 1)
        threshold = threshold * threshold_step
    endwhile
    
    if (peaks_found <= 0)
        return;
    end

    num_peaks = peaks_found;
    if (peaks_found > target_count)
        num_peaks = target_count;
        [sorted_strengths, sorting_perm] = sort(-peaks(:,2));
        permuted_indices = peaks(sorting_perm,1);
        peaks = permuted_indices(1:num_peaks);
    endif

    % +1: Line up peak location indices with spectrum.
    peak_locations = peaks(:,1) + 1;
    % Quadratic interpolation to estimate inter-bin peak locations
    c0 = mag_spectrum(peak_locations);
    to_the_left = mag_spectrum(peak_locations-1);
    to_the_right = mag_spectrum(peak_locations+1);
    c1 = 0.5 * (to_the_right - to_the_left);
    c2 = 0.5 * (to_the_right + to_the_left) - c0;
    xs_at_peaks = -c1 ./ (2 * c2);
    peak_strengths = (c2 .* xs_at_peaks + c1) .* xs_at_peaks + c0;
    % -1: Indices need to be zero-based, not one-based, because the first bin corresponds
    %     to 0 Hz.
    %partials_list = [4096/81*8 1024; 4096/64*8 1024; 4096/49*8 1024];
    partials_list = [xs_at_peaks+peak_locations-1 peak_strengths];
    
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [basis, whole_grain_indices] = build_basis(section, spectrum, partials_list, ...
                                                    window_length, max_grain_length, ...
                                                    time_resolution)

    partial_cycs_per_samp = partials_list(:,1) / size(spectrum,1);
    
    % This creates a 2-D array, columns of complex exponential functions, different
    % frequencies per column
    % Broadcasting is intended here.
    probes = exp(-i * 2 * pi * partial_cycs_per_samp * [0:window_length-1]) .* ...
             (0.5 * (1.0 - cos([0:window_length-1] * 2 * pi / window_length)));
    % Probe length is how much of the array that corresponds to a whole number of cycles
    whole_cycles = floor(partial_cycs_per_samp * window_length);
    probe_lengths = floor(whole_cycles / partial_cycs_per_samp);
    
    responses = probes * section;
    phases = atan2(imag(responses), real(responses))
    
    % For each partial, use the phase to decide where the first crest occurs. This is
    % also the start point of the first grain.
    %  INT ei -wt ei w(t+f) dt
    %  = INT ei wf dt = A ei wf
    %  wf = rads/sec * sec = rads
    %  f = arg/w
    %  ei w(t+f) = 1
    %  w(t+f) = 2 pi k
    %  t = (2 pi k) / w - f
    %  so use smallest k such that t is not negative
    %  2 pi k / w - f >= 0
    %  2 pi k >= wf
    %  k >= wf/2pi = ceil(wf/2pi)
    %  t = (2pi ceil(arg/2pi) - arg) / w
    cycle_fracs = phases / (2*pi);
    first_crests = (ceil(cycle_fracs) - cycle_fracs) ./ partial_cycs_per_samp + 1;
    
    % Figure out wavenumbers
    grain_wavenumbers = floor(partial_cycs_per_samp * max_grain_length / 4) * 4;
    stride_wavenumbers = grain_wavenumbers/2;
    for n = 1:size(stride_wavenumbers,1)
        mod_values = mod(grain_wavenumbers(n), [1:grain_wavenumbers(n)/2]');
        nonzeros = sign(mod_values);
        divisors = find(1-nonzeros);
        max_stride = max(floor(partial_cycs_per_samp * time_resolution));
        not_too_large = (divisors <= max_stride) .* divisors;
        stride_wavenumbers(n) = max(not_too_large)
    endfor

    % For each partial, determine the number of grains that will fit in window_length
    grain_lengths = grain_wavenumbers ./ partial_cycs_per_samp;
    grain_strides = stride_wavenumbers ./ partial_cycs_per_samp;
        
    % x(k) = kL/W + f
    % x(k) + L > 1      AND  x(k) < N
    % kL/W + L + f > 1  AND  kL/W + f < N
    % k > (1-L-f)W/L  AND  k < (N-f)W/L
    % k > floor((1-L-f)W/L) + 1  AND  k < ceil((N-f)W/L) - 1
    % k has (ceil((N-f)W/L) - 1) - (floor((1-L-f)W/L) + 1) + 1 solutions
    % = ceil((N-f)W/L) - 1 - floor((1-L-f)W/L) - 1 + 1
    % = ceil((N-f)W/L) - floor((1-L-f)W/L) - 1
    % = ceil((N-f)W/L) - floor((1-f)W/L - W) - 1
    % = ceil((N-f)W/L) - floor((1-f)W/L) + W - 1
    % = ceil((N-f)W/L) - floor((1-f)W/L)
    % Note that W/L = partial_cycs_per_samp
    % We finally add a division by stride_wavenumbers
    % Things are broken if stride_wavenumbers doesn't divide grain_wavenumbers
    grain_counts = ceil ((window_length-first_crests) .* partial_cycs_per_samp ./ ...
                         stride_wavenumbers) - ...
                   floor((            1-first_crests) .* partial_cycs_per_samp ./ ...
                         stride_wavenumbers) + ...
                   grain_wavenumbers./stride_wavenumbers - 1;
    
    % Determine dimensions of basis matrix
    basis_rows = window_length;
    basis_cols = sum(grain_counts);
    basis = zeros(basis_rows, basis_cols);
    
    % Put stuff in basis matrix
    column_number = 1;
    total_grains = sum(grain_counts);
    whole_grains = zeros(total_grains,1);
    for k = 1:size(partials_list,1)
    
        for j = 1:grain_counts(k)
        
            grain_start = first_crests(k) + grain_strides(k)*(j-1) - grain_lengths(k);
            grain_discrete_start = floor(grain_start) + 1;
            grain_start_frac = grain_discrete_start - grain_start;
            grain_end = grain_start + grain_lengths(k);
            grain_discrete_end = ceil(grain_end) - 1;
            grain_function_start = grain_start_frac;
            grain_function_end = grain_discrete_end - grain_start;
            % Rounding wouldn't be necessary here if we had infinite precision. The
            % value we're rounding is already supposed to be an integer, but, you know,
            % roundoff error. By the way, you can prove it's supposed to be an integer if
            % you do all the substitutions and let the pluses and minuses cancel.
            grain_discrete_length = round(grain_function_end - grain_function_start) + 1;
            grain_function_ordinates = [0:grain_discrete_length]' + grain_function_start;
            grain_raw_function = cos(grain_function_ordinates * 2 * pi * ...
                                     grain_wavenumbers(k) / grain_lengths(k));
            grain_window_function = 0.5 * (1 - cos(grain_function_ordinates * ...
                                                   2 * pi / grain_lengths(k)));
                                
            basis_vector_start = max(grain_discrete_start, 1);
            basis_vector_end   = min(grain_discrete_end, size(basis,1));
            function_index_start = basis_vector_start - grain_discrete_start + 1;
            function_index_end   = basis_vector_end - grain_discrete_start + 1;
            basis(basis_vector_start:basis_vector_end, column_number) = ...
                 grain_raw_function   (function_index_start:function_index_end) .* ...
                 grain_window_function(function_index_start:function_index_end);
            
            if (grain_discrete_start >= 1) && (grain_discrete_end <= size(basis,1))
                whole_grains(column_number) = 1;
            endif
            
            column_number = column_number + 1;
                
        endfor
        
    endfor
    
    whole_grain_indices = find(whole_grains);
    
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
