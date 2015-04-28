
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
    stride_size = window_length/2;
    num_strides = floor((padded_input_size-window_length) / stride_size) + 1;
    for stride = 0:num_strides-1
        % Section goes from here to here + window length
        section_begin = stride * stride_size + 1;
        section_end = section_begin + window_length - 1;
        section = padded_input(section_begin:section_end);
        
        residual = section;
        
        oversample_exp = 3;
        
        % Do the following until there don't seem to be any salient partials
        %while (0 == 0)

            % Take windowed FFT
            window_function = lorentz(window_length);
            % window_function = 0.5 * (1 + cos([-window_length/2+1:window_length/2]' * ...
            %                                 2 * pi / window_length));
            windowed_section = residual .* window_function;
            spectrum = superfft(windowed_section, oversample_exp);

            % Find salient partials
            partials_list = find_partials(spectrum, 63);
            partials_list(:,1) = partials_list(:,1) / 2^oversample_exp;
keyboard();
            
            % Build basis vectors from all partials
            max_grain_length = window_length / 10;
            time_resolution = 100;
            [basis whole_grain_indices] = ...
                build_basis(residual, partials_list, ...
                            max_grain_length, time_resolution);
keyboard();
            
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
        
        padded_input(section_begin:section_end) = residual;
        
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

function partials_list = find_partials(spectrum, mwidth)
    
    spectrum = spectrum(:);
    
    % Enforce odd mwidth
    mwidth = floor(mwidth/2) * 2 + 1;

    mag_spectrum = abs(spectrum);
    half_spectrum = mag_spectrum(1 : end/2);
    half_spectrum_size = size(half_spectrum, 1);
    
    moving_window = zeros(half_spectrum_size-mwidth+1, mwidth);
    for k = 1:mwidth
        moving_window(:,k) = half_spectrum(k:end-mwidth+k, 1);
    endfor
    
    med_filt_output = median(moving_window, 2);
    % Many values in salient_values will be zero due to the subtraction here.
    salient_values_start = (mwidth+1)/2;
    salient_values_end = -salient_values_start + 1;
    salient_values = half_spectrum(salient_values_start:end+salient_values_end) - ...
                     med_filt_output;
    % Because of all the zeros in salient_values, we need to deliberately reject peaks
    % at height zero, hence the third factor in this multiplication.
    % The fourth and fifth factors check that these indices truly correspond to peaks in
    % the original spectrum. Sometimes a little lump on a hill can fool the algorithm
    % otherwise.
    indicators_of_pos_peaks = ...
        (salient_values(3:end) <= salient_values(2:end-1)) .* ...
        (salient_values(2:end-1) > salient_values(1:end-2)) .* ...
        (salient_values(2:end-1) > 0) .* ...
        (half_spectrum(salient_values_start+2:end+salient_values_end) <= ...
         half_spectrum(salient_values_start+1:end+salient_values_end-1)) .* ...
        (half_spectrum(salient_values_start+1:end+salient_values_end-1) > ...
         half_spectrum(salient_values_start:end+salient_values_end-2));
    peak_indices = find(indicators_of_pos_peaks);
    %[sorted_peak_heights peak_perm] = sort(salient_values(peak_indices));
    %accum_energy_by_freq = cumsum(sqrt(sorted_peak_heights .* sorted_peak_heights));
    %% We're trying to find a few peaks that account for more than their share of the energy
    %% in the spectrum. normalized_accum is a non-decreasing series with a non-decreasing
    %% derivative that captures how flat the distribution is. (It's something like the
    %% integral of a histogram.) A purely flat spectrum will result in a flat diagonal. The
    %% more "peaky" the distribution, the more prominent an elbow there will be. Our ad hoc
    %% criteria are to say that the elbow is the point where the slope exceeds 1.
    %% It probably is possible to do this with the derivative of this function, but we're not
    %% going to put the time in for that right now, as this approach already works.
    %crook_of_elbow = ...
       %find((accum_energy_by_freq(3:end) - accum_energy_by_freq(1:end-2)) > 2)(1);

    %% The sorting permutation tells us which elements in peak_indices these peaks came from.
    %strongest_peak_indices = peak_perm(crook_of_elbow:end);
    %% But peak_indices is a subset of salient_values. So we need to go backward one more
    %% step and figure out which bins in salient_values correspond to these elements in
    %% peak_indices.
    %strongest_peak_bins = peak_indices(peak_perm);
    % Adjust these indices, which point into salient_values, so they point properly into
    % half_spectrum.
    %strongest_peak_bins = strongest_peak_bins + salient_values_start;
    strongest_peak_bins = peak_indices + salient_values_start;
    c0 = half_spectrum(strongest_peak_bins);
    to_the_left = half_spectrum(strongest_peak_bins-1);
    to_the_right = half_spectrum(strongest_peak_bins+1);
    c1 = 0.5 * (to_the_right - to_the_left);
    c2 = 0.5 * (to_the_right + to_the_left) - c0;
    xs_at_peaks = -c1 ./ (2 * c2);
    peak_strengths = (c2 .* xs_at_peaks + c1) .* xs_at_peaks + c0;
    % -1: Indices need to be zero-based, not one-based, because the first bin corresponds
    %     to 0 Hz.
    partials_list = [xs_at_peaks+strongest_peak_bins-1 peak_strengths];
keyboard();
    
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [basis, whole_grain_indices] = build_basis(section, partials_list, ...
                                                    max_grain_length, time_resolution)

    section_length = size(section,1);
    partial_cycs_per_samp = partials_list(:,1) / section_length;
    
    % This creates a 2-D array, columns of complex exponential functions, different
    % frequencies per column
    % Broadcasting is intended here.
    probes = exp(-i * 2 * pi * partial_cycs_per_samp * [0:section_length-1]) .* ...
             (0.5 * (1.0 - cos([0:section_length-1] * 2 * pi / section_length)));
    % Probe length is how much of the array that corresponds to a whole number of cycles
    whole_cycles = floor(partial_cycs_per_samp * section_length);
    probe_lengths = floor(whole_cycles / partial_cycs_per_samp);
    
    responses = probes * section;
    phases = atan2(imag(responses), real(responses));
    
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
keyboard();
    grain_wavenumbers = floor(partial_cycs_per_samp * max_grain_length / 4) * 4;
    grain_wavenumbers = grain_wavenumbers(find(grain_wavenumbers));
    if (size(grain_wavenumbers,1) == 0)
        basis = [];
        whole_grain_indices = [];
        return;
    endif
    stride_wavenumbers = grain_wavenumbers/2;
    for n = 1:size(stride_wavenumbers,1)
        mod_values = mod(grain_wavenumbers(n), [1:grain_wavenumbers(n)/2]');
        nonzeros = sign(mod_values);
        divisors = find(1-nonzeros);
        max_stride = max(floor(partial_cycs_per_samp * time_resolution));
        not_too_large = (divisors <= max_stride) .* divisors;
        stride_wavenumbers(n) = max(not_too_large);
    endfor

    % For each partial, determine the number of grains that will fit in section_length
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
    grain_counts = ceil ((section_length-first_crests) .* partial_cycs_per_samp ./ ...
                         stride_wavenumbers) - ...
                   floor((            1-first_crests) .* partial_cycs_per_samp ./ ...
                         stride_wavenumbers) + ...
                   grain_wavenumbers./stride_wavenumbers - 1;
    
    % Determine dimensions of basis matrix
    basis_rows = section_length;
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
