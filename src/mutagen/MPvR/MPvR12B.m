
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grains, residual] = MPvR(in, window_length, oversample_exp, noise_floor, depth)

    grains = zeros(0);
    
    in = in(:);
    in_length = size(in, 1);

    % Zero-pad input so its length is a multiple of window_length/2

    step_size = window_length / 2;

    % We use step_size*2 every step, so we stop at the penultimate step.
    % That's why we subtract 1 here.
    num_steps = ceil(in_length / step_size) - 1;
    
    % We use step_size*2 every step, so we stop at the penultimate step.
    % That's why we subtract 1 here.
    num_steps = ceil(in_length / step_size) - 1;

    padded_length = (num_steps+1) * step_size;
    corrected_input = [in ; zeros(padded_length-in_length, 1)];
    oversample_factor = 2^oversample_exp;

    % For each section
    for loop_number = 0:num_steps-1

        % Grab window's worth
        section_begin = 1 + loop_number * step_size;
        section_end = section_begin + window_length - 1;
        section = corrected_input(section_begin:section_end);

        % Apply window.
        window_function = 0.5 * (1 + cos([-window_length/2+1:window_length/2]' * ...
                                         2 * pi / window_length));
        windowed_section = section .* window_function;

        while (0 == 0)

            % Find spectral peaks
            complex_spectrum = superfft(windowed_section, oversample_exp);
            peak_couples = find_peaks(complex_spectrum, 1)
            
            num_peaks = size(peak_couples, 1);
            if (num_peaks > 0)
                % Correct the indices to account for frequency oversampling.
                peak_couples(:,1) = (peak_couples(:,1) - 1) / oversample_factor + 1;
            endif

            % Remove energy at each peak. We could also do just the strongest peak, then
            % go back and re-evaluate the spectrum now that it has been altered, in which
            % case we wouldn't have a for loop here, just one pass that handles the strongest
            % peak returned by the peak finder.
            for peak_number = 1:num_peaks
            
                % We still need to subtract 1 to get the proper frequency because
                % the indexing of all these vectors starts at 1, but the spectrum
                % starts at 0 Hz.
                peak_freq_rads_per_samp = ...
                    2 * pi * (peak_couples(peak_number,1) - 1) / window_length;

                % Decompose convolution of kernel and section into a train of kernels
                [new_grains residual] = ...
                    extract_frequency(section, peak_freq_rads_per_samp, 1e300);
                section = residual;
                
                % Add kernels to list
                % The new grains come to us in a time coordinate system that's different from
                % the one we're using here. So they need to be translated.
                new_grains(:,1) = new_grains(:,1) + section_begin;
                grains = [grains ; new_grains];
                
            endfor
            
            % Perhaps go back and re-evaluate this section until we reach a threshold?
            % Or just break out. If we hard-code breaking out, there's no need for the
            % encasing while loop
            break;
            
        endwhile

        corrected_input(section_begin:section_end) = section;
        
        % We've found all the partials to be found at this scale, so recurse.
        % This will complete the decomposition of this section.

        % Step to the next window. Remember that window overlap is 50%.
        
    endfor

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = upsample(in, upsample_exp)

    in = in(:);
    
    if (upsample_exp < 1)
        out = in;
        return;
    endif
    
    in_length = size(in,1);
    for e = 1:31
        if (2^e > in_length)
            padded_in_length = 2^e;
            break;
        endif
    endfor
    
    padded_in = [in ; zeros(padded_in_length-in_length,1)];
    spectrum = fft(padded_in);
    upsample_factor = 2^upsample_exp;
    spectrum_of_upsample = [spectrum(1:size(spectrum,1)/2-1) ; ...
                            zeros(size(spectrum,1) * (upsample_factor-1), 1) ; ...
                            spectrum(size(spectrum,1)/2:size(spectrum,1))];

    upsampled_padded_in = ifft(spectrum_of_upsample);
    out = upsampled_padded_in(1 : in_length*upsample_factor);

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

function [c3,c2,c1,c0] = cubic_interp(xm1,xp0,xp1,xp2)
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function array_of_peaks = find_peaks(spectrum, min_deviation)

    array_of_peaks = zeros(0);
    
    % Compute magnitudes
    spectrum = spectrum(:);
    mag_spectrum_length = size(spectrum, 1) / 2 + 1;
    mag_spectrum = abs(spectrum(1:mag_spectrum_length));

    % Take the first difference.
    first_diff = [mag_spectrum(1) - mag_spectrum(2) ; ...
                  mag_spectrum(2:end) - mag_spectrum(1:end-1) ; ...
                  mag_spectrum(end-1) - mag_spectrum(end)];
    signs = sign(first_diff);
    % Some of these signs might be zero. The usual peak's first difference looks like +-,
    % but we will also accept the +0 sequence. So convert all zeros to negatives.
    signs = signs - 0.5;
    signs = sign(signs);
    % We can use this info to get zero crossings of the first difference.
    peakiness = signs(1:end-1) - signs(2:end);
    % Each value will be either -2, 0, or 2. After this the possible values will be 0 or 8.
    % The nonzeros will be the locations of the peaks.
    non_peaks_removed = peakiness .* (peakiness + 2);
    % This will find the locations of peaks.
    peak_locations = find(non_peaks_removed);

    if (size(peak_locations, 1) <= 0)
        return;
    end

    standard_deviation = std(mag_spectrum);
    peak_strengths = mag_spectrum(peak_locations) / standard_deviation;
    above_threshold = sign(sign(peak_strengths - min_deviation) - 0.5);
    above_threshold = above_threshold .* (above_threshold + 1);
    salient_peak_locations = peak_locations(find(above_threshold));
    
    if (size(salient_peak_locations, 1) <= 0)
        return;
    end
    
    extended_mag_spectrum = [mag_spectrum(2) ; mag_spectrum ; mag_spectrum(end-1)];
    % Quadratic interpolation to estimate inter-bin peak locations
    c0 = extended_mag_spectrum(salient_peak_locations+1);
    to_the_left = extended_mag_spectrum(salient_peak_locations);
    to_the_right = extended_mag_spectrum(salient_peak_locations + 2);
    c1 = 0.5 * (to_the_right - to_the_left);
    c2 = 0.5 * (to_the_right + to_the_left) - c0;
    x_at_peak = -c1 ./ (2 * c2);
    peak = (c2 .* x_at_peak + c1) .* x_at_peak + c0;
    array_of_peaks = [x_at_peak+salient_peak_locations peak];
    
    % Sort the peak values column, then permute the peak indices column
    % accordingly.
    [array_of_peaks(:,2), perm] = sort(array_of_peaks(:,2), "descend");
    array_of_peaks(:,1) = array_of_peaks(perm,1);

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grains, residual] = ...
    extract_frequency(section, rads_per_sample, threshold_significance)

    grains = zeros(0,4);
    
    analysis_wavenum = 16;
    decomp_wavenum = 4;

    % Build analysis kernel
    % kernel_start and kernel_end are integers, opposite signs, which means that
    % the length of the kernel will be odd and the midpoint will be at 0.
    analysis_kernel_start = 1 - ceil(analysis_wavenum * pi / rads_per_sample);
    analysis_kernel_end = -analysis_kernel_start;
    analysis_kernel_sinus = cos([analysis_kernel_start:analysis_kernel_end]' * rads_per_sample);
    analysis_kernel_window = ...
        0.5 * (1 + cos([analysis_kernel_start:analysis_kernel_end]' * rads_per_sample / ...
                       analysis_wavenum));
    analysis_kernel = analysis_kernel_sinus .* analysis_kernel_window;
    analysis_kernel_amplitude = 1 / sqrt(sum(analysis_kernel .* analysis_kernel));
    analysis_kernel = analysis_kernel * analysis_kernel_amplitude;
    analysis_kernel_discrete_length = size(analysis_kernel, 1); % Must be odd

    % Convolve
    filter_output = fftconv(section, analysis_kernel);
    filter_output_length = size(filter_output, 1);
    
    % Build reference decomposition kernel
    % Length of this kernel will be odd. See analysis kernel.
    ref_decomp_kernel_start = 1 - ceil(decomp_wavenum * pi / rads_per_sample);
    ref_decomp_kernel_end = -ref_decomp_kernel_start;
    ref_decomp_kernel_sinus = cos([ref_decomp_kernel_start:ref_decomp_kernel_end]' * ...
                                  rads_per_sample);
    ref_decomp_kernel_window = ...
        0.5 * (1 + cos([ref_decomp_kernel_start:ref_decomp_kernel_end]' * rads_per_sample / ...
                       decomp_wavenum));
    ref_decomp_kernel = ref_decomp_kernel_window .* ref_decomp_kernel_sinus;
    ref_decomp_kernel_amplitude = 1 / sqrt(sum(ref_decomp_kernel .* ref_decomp_kernel));
    ref_decomp_kernel = ref_decomp_kernel * ref_decomp_kernel_amplitude;
    ref_decomp_kernel_discrete_length = size(ref_decomp_kernel, 1);

    % Build analysis decomposition kernel
    ad_kernel = fftconv(analysis_kernel, ref_decomp_kernel);
    ad_kernel = ad_kernel / max(ad_kernel);
    ad_kernel_discrete_length = size(ad_kernel, 1);

    salient_filt_oput_begin = analysis_kernel_discrete_length + ...
                              (ref_decomp_kernel_discrete_length-1)/2;
    last_peak_index = salient_filt_oput_begin;
    salient_filt_oput_end = size(section,1) - analysis_kernel_discrete_length - ...
                            (ref_decomp_kernel_discrete_length-1)/2;

    while 0 == 0
        % Find next peak, or skip out if there is none
        found_peak = 0;
        last_peak_int_index = salient_filt_oput_begin;
        for i = last_peak_int_index+1:salient_filt_oput_end
            if (filter_output(i-1) < filter_output(i) && filter_output(i) >= filter_output(i+1))
                last_peak_int_index = i;
                c0 = filter_output(i);
                to_the_left = filter_output(i-1);
                to_the_right = filter_output(i+1);
                c1 = 0.5 * (to_the_right - to_the_left);
                c2 = 0.5 * (to_the_right + to_the_left) - c0;
                last_peak_frac_index = -c1 ./ (2 * c2);
                last_peak = (c2 * last_peak_frac_index + c1) * last_peak_frac_index + c0;
                found_peak = 1;
                break;
            endif
        endfor
        
        if (found_peak != 1)
            break;
        endif
        
        % Find maximum amplitude scale of AD kernel that won't screw things up
        
        %  Locate peaks of absolute AD kernel
        ad_peak_locations = [0:pi/rads_per_sample:ad_kernel_discrete_length]' + 1;
        ad_peak_locations = ad_peak_locations(2:end-1);
        ad_peak_int_locations = round(ad_peak_locations);
        ad_peak_frac_locations = ad_peak_locations - ad_peak_int_locations;
        c0 = ad_kernel(ad_peak_int_locations);
        to_the_left = ad_kernel(ad_peak_int_locations-1);
        to_the_right = ad_kernel(ad_peak_int_locations+1);
        c1 = 0.5 * (to_the_right - to_the_left);
        c2 = 0.5 * (to_the_right + to_the_left) - c0;
        ad_peaks = (c2 .* ad_peak_frac_locations + c1) .* ad_peak_frac_locations + c0;
        
        %  Find the peak in the AD kernel whose index is closest to
        %  (decomp_kernel_discrete_length-1)/2
        ad_peak_locations = ad_peak_int_locations + ad_peak_frac_locations;
        [alignment_ad_peak, alignment_ad_peak_index] = ...
            min(abs(ad_peak_locations - (ref_decomp_kernel_discrete_length-1)/2));
keyboard();        
        %  For each AD kernel peak, determine the amplitude of the AD kernel that will send
        %   the corresponding filter output value to 0
        last_peak_index = last_peak_int_index + last_peak_frac_index;
        scale_analysis_locations = ...
            ad_peak_locations - ad_peak_locations(alignment_ad_peak_index) + last_peak_index;
        scale_analysis_int_locations = round(scale_analysis_locations);
        c0 = filter_output(scale_analysis_int_locations);
        to_the_left = filter_output(scale_analysis_int_locations-1);
        to_the_right = filter_output(scale_analysis_int_locations+1);
        c1 = 0.5 * (to_the_right - to_the_left);
        c2 = 0.5 * (to_the_right + toad_the_left) - c0;
        scale_analysis_frac_locations = scale_analysis_locations - scale_analysis_int_locations;
        filt_values_at_ad_peaks = ...
            (c2 .* scale_analysis_frac_locations + c1) .* scale_analysis_frac_locations + c0;
                
        %  Select the least such amplitude
        min_ad_scale = max(min(filt_values_at_ad_peaks ./ ad_peaks), 0);
        
        % Subtract ad_kernel from filter_output
        interpolating_frac_index = ad_peak_frac_locations(alignment_ad_peak_index) - ...
                                   last_peak_frac_index;
        c0 = ad_kernel;
        to_the_left = [0 ; ad_kernel(1:end-1)];
        to_the_right = [ad_kernel(2:end) ; 0];
        c1 = 0.5 * (to_the_right - to_the_left);
        c2 = 0.5 * (to_the_right + to_the_left) - c0;
        interpolated_ad_kernel = ...
            (c2 .* interpolating_frac_index + c1) .* interpolating_frac_index + c0;
        ad_removal_start = ceil(last_peak_index - (decomp_kernel_discrete_length-1)/2);
        ad_removal_end = ad_removal_start + ad_kernel_discrete_length - 1;
        filter_output(ad_removal_start:ad_removal_end) = ...
            filter_output(ad_removal_start:ad_removal_end) - ...
            min_ad_scale * interpolated_ad_kernel;
        
        % Construct decomposition kernel
        decomp_kernel_start_time = 0;
        if (interpolating_frac_index > 0)
            decomp_kernel_start = -(decomp_kernel_discrete_length-1)/2 - ...
                                  interpolating_frac_index;
        else
            decomp_kernel_start = -(decomp_kernel_discrete_length-1)/2 - ...
                                  interpolating_frac_index + 1;
        endif
        decomp_kernel_end_time = decomp_kernel_start + decomp_kernel_discrete_length;
        if (interpolating_frac_index == 0)
            decomp_kernel_end = decomp_kernel_end - 1;
        endif
        decomp_kernel_sinus = cos([decomp_kernel_start:decomp_kernel_end]' * ...
                                  rads_per_sample);
        decomp_kernel_window = ...
           0.5 * (1 + cos([decomp_kernel_start:decomp_kernel_end]' * rads_per_sample / ...
                          decomp_wavenum));
        decomp_kernel_amplitude = 1 / sqrt(sum(ref_decomp_kernel .* ref_decomp_kernel));
        decomp_kernel = decomp_kernel * decomp_kernel_amplitude * min_ad_scale;
        decomp_kernel_discrete_length = size(decomp_kernel,1);
        
        % Subtract decomp_kernel from section
        removal_start = ad_removal_start + analysis_kernel_discrete_length - 1;
        removal_end = removal_start + decomp_kernel_discrete_length;
        section(removal_start:removal_end) = section(removal_start:removal_end) - ...
                                             decomp_kernel;
        
    endwhile

    residual = section;

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
