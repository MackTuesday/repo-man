
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grains, residual] = MPvR(in, window_length, oversample_exp, noise_floor, depth)

    grains = 0;
    residual = zeros(window_length*2, 1);

    % Zero-pad input so its length is a multiple of window_length/2

    step_size = window_length / 2;

    % We use step_size*2 every step, so we stop at the penultimate step.
    % That's why we subtract 1 here.
    num_steps = ceil(in_length / step_size) - 1;

    padded_length = (num_steps+1) * step_size;
    corrected_input = [in ; zeros(padded_length-in_length, 1)];
    oversample_factor = 2^oversample_exp;

    for loop_number = 0:num_steps-1

        % Grab window's worth
        section_begin = 1 + loop_number * step_size;
        section_end = section_begin + window_length - 1;
        section = corrected_input(section_begin:section_end);
if (depth<2)
disp([depth section_begin section_end]);
endif
        % Apply window. The first half is ones because we assume the previous
        % window, which tapers off in the second half, made it so the section
        % we're working with now should already be tapered up in the first half,
        % because we supposedly decomposed it and now it's totally subtracted
        % away from the original signal. If it isn't tapered up, we're at the
        % very first part of the whole signal. A lack of tapering at the
        % beginning constitutes a transient, which we don't want to miss anyway.
        if (loop_number == num_steps-1)
            window_function = ones(window_length,1);
        else
            window_function = 0.5 * (1 + cos([1:window_length/2]' *
                                             2 * pi / window_length));
            window_function = cat(1, ones(window_length/2,1), window_function);
        endif
        windowed_section = section .* window_function;

        while (0 == 0)
        
            % Get a list of valid peaks in the spectrum, sorted by strength.
            complex_spectrum = superfft(windowed_section, oversample_exp);
            peak_couples = find_peaks(complex_spectrum, 1);
if depth <= 3
figure(1)
plot(windowed_section);
figure(2)
plot(abs(complex_spectrum));
peak_couples
keyboard();
endif

            num_peaks = size(peak_couples, 1);
            if (num_peaks > 0)
                % Correct the indices to account for frequency oversampling.
                peak_couples(:,1) = (peak_couples(:,1) - 1) / oversample_factor + 1;
            endif

            % Loop through them, trying to remove one.
            removed_peak = 0;
            peak_number = 1;
            for peak_number = 1:num_peaks
                % Build analysis kernel. We make it half the section length
                % because the window function on the kernel and the window
                % function on the analysis section (constructed below) must
                % overlap 100%. If the window functions are the same size, we get
                % no time resolution because the window functions would coincide
                % exactly.

                % We still need to subtract 1 to get the proper frequency because
                % the indexing of all these vectors starts at 1, but the spectrum
                % starts at 0 Hz.
                peak_freq_rads_per_samp = ...
                    2 * pi * (peak_couples(peak_number,1) - 1) / window_length;
                kernel_length = window_length / 2;
                kernel_window_function = ...
                    0.5 * (1 + cos([-kernel_length/2+1:kernel_length/2]' * ...
                                   2 * pi / kernel_length));
                filter_kernel = cos([-kernel_length/2+1:kernel_length/2]' * ...
                                    peak_freq_rads_per_samp) .* ...
                                kernel_window_function;
                filter_kernel_amplitude = ...
                    1 / sqrt(sum(filter_kernel .* filter_kernel));
                filter_kernel = filter_kernel * filter_kernel_amplitude;

                % Run section through filter.
                filter_output = fftconv(windowed_section, filter_kernel);

                % Find highest absolute peak in the response
                filter_output_length = size(filter_output, 1);
                max_response = 0;
                x_at_max_response = 0;
                frac_x_at_max_response = 0;
                [max_element, int_x_at_max_response] = max(abs(filter_output));
                if (int_x_at_max_response > 1 && ...
                    int_x_at_max_response < filter_output_length)
                    c0 = filter_output(int_x_at_max_response);
                    c1 = 0.5 * (filter_output(int_x_at_max_response+1) - 
                                filter_output(int_x_at_max_response-1));
                    c2 = 0.5 * (filter_output(int_x_at_max_response+1) + 
                                filter_output(int_x_at_max_response-1)) - c0;
                    frac_x_at_max_response = -c1 / (2 * c2);
                    max_response = (c2 * frac_x_at_max_response + c1) * ...
                                   frac_x_at_max_response + c0;
                    x_at_max_response = int_x_at_max_response + ...
                                        frac_x_at_max_response;
                endif
                
                % If the max response is in such a place that the corresponding
                % wavelet doesn't overlap with the windowed section 100%, this
                % peak is no good and we need to keep looking.
                salient_filt_oput_begin = kernel_length;
                salient_filt_oput_end = ...
                    salient_filt_oput_begin + window_length - kernel_length;
                if (x_at_max_response < salient_filt_oput_begin || ...
                    x_at_max_response > salient_filt_oput_end)
                   peak_number = peak_number + 1;
                   continue;
                endif
                
                % This peak is good. We will remove it.
                
                % Translate ordinates so they're based on the windowed
                % section rather than the filter output.
                int_x_max_resp_trans = int_x_at_max_response - kernel_length + 1;
                x_at_max_resp_trans = x_at_max_response - kernel_length + 1;

                % Generate wavelet based on magnitude peak.
                
                % The wavelet has the same shape as the filter kernel. We seek
                % such a scale for it that its square-sum matches the maximum
                % response value we have found. The kernel was constructed with
                % a norm of one, so we multiply the kernel's amplitude by the
                % response to get the norm of the wavelet where we want it.
                wavelet_amplitude = filter_kernel_amplitude * max_response;
                wavelet_length = kernel_length;
                
                % There's no point in letting the start magnitude of the wavelet
                % be 0.
                wavelet_t_first = -frac_x_at_max_response - wavelet_length/2 + ...
                                  ceil(frac_x_at_max_response);
                wavelet_t_last = wavelet_t_first + wavelet_length - 1;
                wavelet_window_function = 0.5 * (1 + cos([wavelet_t_first:wavelet_t_last]' * ...
                                                         2 * pi / wavelet_length));
                wavelet = wavelet_amplitude * ...
                          cos([wavelet_t_first:wavelet_t_last]' * ...
                              peak_freq_rads_per_samp) .* wavelet_window_function;
                % Subtract wavelet from section
                removal_start = int_x_max_resp_trans + ceil(frac_x_at_max_response);
                removal_end = removal_start + wavelet_length - 1;
if depth <= 3
figure(1)
plot(windowed_section(removal_start:removal_end));
endif
                windowed_section(removal_start:removal_end) = ...
                    windowed_section(removal_start:removal_end) - wavelet;
if depth <= 3
figure(2)
plot(wavelet)
figure(3)
plot(windowed_section(removal_start:removal_end));
endif
                % Add wavelet to result list
                new_grain = [x_at_max_resp_trans+section_begin-1 ...
                             peak_freq_rads_per_samp wavelet_amplitude wavelet_length];
if depth <= 3
new_grain
keyboard();
endif
                if (grains == 0)
                    grains = new_grain;
                else
                    grains(end+1,:) = new_grain;
                endif

                % We have removed the peak. Now see about the next one in the
                % list.
                peak_number = peak_number + 1;
                
            end
            
            % We've found all the partials to be found at this scale, so recurse.
            % This will complete the decomposition of this windowed section.
            if (window_length >= 32)
                new_grains = 0;
                [new_grains, new_residual] = ...
                    MPvR(windowed_section, window_length/2, oversample_exp, ...
                         noise_floor, depth+1);

                if (new_grains(1) ~= 0)
                    % The positions of the grains we got back are in the
                    % coordinate system used by the recursed call to MPvR.
                    % They need to be transformed to our coordinate system.
                    % The frequencies are in radians per sample, which is a
                    % global coordinate system. So is the amplitude scale.
                    new_grains(:,1) = new_grains(:,1) + section_begin - 1;
                    
                    if (grains == 0)
                        grains = new_grains;
                    else
                        grains = [grains ; new_grains];
                    endif
                endif
                
                residual = [residual ; zeros(section_end-size(residual,1),1)];
                residual(section_begin:section_end) = ...
                    residual(section_begin:section_end) + new_residual;
            else
                residual = [residual ; zeros(section_end-size(residual,1),1)];
                residual(section_begin:section_end) = ...
                    residual(section_begin:section_end) + windowed_section;
            endif
            
            % The decomposition of this windowed section is complete, so subtract
            % it from the input and move on.
            corrected_input(section_begin:section_end) = ...
                corrected_input(section_begin:section_end) - section .* window_function;
                
            break;
            
        endwhile

        % We have exited from the loop that decomposes the windowed section. We
        % step to the next window. Remember that window overlap is 50%.
        
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
