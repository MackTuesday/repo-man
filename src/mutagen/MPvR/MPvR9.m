
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
            peak_couples = find_peaks(complex_spectrum, 1);
            
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
                new_grains = extract_frequency(section, peak_freq_rads_per_samp);
                
                % Add kernels to list
                grains = [grains ; new_grains];
                
            endfor
            
            % Perhaps go back and re-evaluate this section until we reach a threshold?
            % Or just break out. If we hard-code breaking out, there's no need for the
            % encasing while loop
            break;
            
        endwhile

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
    c1 = 0.5 * (to_the_left - to_the_right);
    c2 = 0.5 * (to_the_left + to_the_right) - c0;
    x_at_peak = -c1 ./ (2 * c2);
    peak = (c2 .* x_at_peak + c1) .* x_at_peak + c0;
    array_of_peaks = [x_at_peak+salient_peak_locations peak];
    
    % Sort the peak values column, then permute the peak indices column
    % accordingly.
    [array_of_peaks(:,2), perm] = sort(array_of_peaks(:,2), "descend");
    array_of_peaks(:,1) = array_of_peaks(perm,1);

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grains, residual] = extract_frequency(section, rads_per_sample)

    grain = zeros(0);
    
    wavenumber = 4;
    
    % Build kernel
    kernel_length = ceil(wavenumber * 2 * pi / rads_per_sample);
    domain_start = -ceil(kernel_length/2) + 1;
    domain_end = domain_start + kernel_length - 1;
    grain_window_function = ...
        0.5 * (1 + cos([domain_start:domain_end]' * 2 * pi / kernel_length));
    filter_window_function = grain_window_function .* grain_window_function;
    filter_kernel_sinus = cos([domain_start:domain_end]' * rads_per_sample);
    filter_kernel = filter_kernel_sinus .* grain_window_function;
    filter_kernel_amplitude = 1 / sqrt(sum(filter_kernel .* filter_kernel));
    % Yes, this is weird. It has to do with matching up grain amplitude with filter responses.
    filter_kernel = filter_kernel .* grain_window_function * filter_kernel_amplitude;
    
    % Convolve
    filter_output = fftconv(section, filter_kernel);
    filter_output_length = size(filter_output, 1);
    
    % While the energy of the residual is greater than some threshold
    residual = section;
    %salient_filt_oput_begin = kernel_length;
    %salient_filt_oput_end = salient_filt_oput_begin + size(section,1) - kernel_length;
    while (0 == 0)
    
        % Identify peak response.
        output_diffs = filter_output(2:end) - filter_output(1:end-1);
        signs_of_output_diffs = sign(sign(output_diffs) - 0.5);
        zero_crossings = signs_of_output_diffs(2:end) - signs_of_output_diffs(1:end-1);
        peak_locations = find(zero_crossings) + 1;
        [max_element, nth_peak] = max(filter_output(peak_locations));
        int_x_at_max_response = peak_locations(nth_peak);
        c0 = filter_output(int_x_at_max_response);
        c1 = 0.5 * (filter_output(int_x_at_max_response+1) - ...
                    filter_output(int_x_at_max_response-1));
        c2 = 0.5 * (filter_output(int_x_at_max_response+1) + ...
                    filter_output(int_x_at_max_response-1)) - c0;
        frac_x_at_max_response = -c1 / (2 * c2);
        max_response = (c2 * frac_x_at_max_response + c1) * frac_x_at_max_response + c0;
        x_at_max_response = int_x_at_max_response + frac_x_at_max_response;
        
        % The grain has the same shape as the filter kernel. We seek
        % such a scale for it that its square-sum matches the maximum
        % response value we have found. The kernel was constructed with
        % a norm of one, so we multiply the kernel's amplitude by the
        % response to get the norm of the wavelet where we want it.
        % You might notice that the grain uses a different envelope from
        % the filter kernel and so maybe this is a mismatch. It isn't.
        % Hint: The filter kernel envelope is deliberately chosen to be
        % the square of the grain envelope.
        grain_amplitude = filter_kernel_amplitude * max_response;
        grain_length = kernel_length;  
                
        % Build grain
        grain_t_first = -frac_x_at_max_response - grain_length/2 + ...
                        ceil(frac_x_at_max_response);
        grain_t_last = grain_t_first + grain_length - 1;
        grain = grain_amplitude * cos([grain_t_first:grain_t_last]' * rads_per_sample) .* ...
                grain_window_function;
                              
        % Translate ordinates so they're based on the windowed
        % section rather than the filter output.
        int_x_max_resp_trans = int_x_at_max_response - kernel_length + 1;
        x_at_max_resp_trans = x_at_max_response - kernel_length + 1;

        % Subtract grain from residual
        removal_start = int_x_max_resp_trans + ceil(frac_x_at_max_response);
        removal_start_correction = 0;
        if (removal_start < 1)
            removal_start_correction = 1 - removal_start;
        endif
        removal_end = removal_start + grain_length - 1;
        removal_end_correction = 0;
        if (removal_end > size(residual,1))
            removal_end_correction = size(residual,1) - removal_end;
        endif
        
        section(removal_start+removal_start_correction: ...
                removal_end+removal_end_correction) = ...
            section(removal_start+removal_start_correction: ...
                    removal_end+removal_end_correction) - ...
            grain(1+removal_start_correction:end+removal_end_correction);
figure(1);
plot(section);
%figure(2);
%plot(grain(1+removal_start_correction:end+removal_end_correction));
keyboard();
        
        % Add grain to list
        if (removal_start_correction == 0 && removal_end_correction == 0)
            residual(removal_start:removal_end) = ...
                residual(removal_start:removal_end) - grain;
            new_grain = [x_at_max_resp_trans-1 rads_per_sample grain_amplitude grain_length];
            grains(end+1,:) = new_grain;
        endif
        
        % Convolve again
        filter_output = fftconv(section, filter_kernel);

    endwhile

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
