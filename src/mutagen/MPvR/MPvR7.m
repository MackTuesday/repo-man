
###############################################################################################

function [grains, residual] = MPvR(in, window_length, oversample_exp, noise_floor, depth)

    grains = 0;
    residual = zeros(window_length*2, 1);

    in = in(:);
    in_length = size(in, 1);

    # Zero-pad input so its length is a multiple of window_length/2

    step_size = window_length / 2;

    # We use step_size*2 every step, so we stop at the penultimate step.
    # That's why we subtract 1 here.
    num_steps = ceil(in_length / step_size) - 1;

    padded_length = (num_steps+1) * step_size;
    corrected_input = [in ; zeros(padded_length-in_length, 1)];
    oversample_factor = 2^oversample_exp;

    for loop_number = 0:num_steps-1

        # Grab window's worth
        section_begin = 1 + loop_number * step_size;
        section_end = section_begin + window_length - 1;
        section = corrected_input(section_begin:section_end);
if (depth<2)
disp([depth section_begin section_end]);
endif
        # Apply window. The first half is ones because we assume the previous
        # window, which tapers off in the second half, made it so the section
        # we're working with now should already be tapered up in the first half,
        # because we supposedly decomposed it and now it's totally subtracted
        # away from the original signal. If it isn't tapered up, we're at the
        # very first part of the whole signal. A lack of tapering at the
        # beginning constitutes a transient, which we don't want to miss anyway.
        if (loop_number == num_steps-1)
            window_function = ones(window_length,1);
        else
            window_function = 0.5 * (1 + cos([1:window_length/2]' *
                                             2 * pi / window_length));
            window_function = cat(1, ones(window_length/2,1), window_function);
        endif
        windowed_section = section .* window_function;

        while (0 == 0)
        
            # Get a list of valid peaks in the spectrum, sorted by strength.
            complex_spectrum = superfft(windowed_section, oversample_exp);
keyboard();
            peak_couples = find_peaks(complex_spectrum, noise_floor, ...
                                      ceil(2.5 * oversample_factor), ...
                                      8 * oversample_factor + 1);

            num_peaks = size(peak_couples, 1);
            if (peak_couples(1,1) == 0)
                num_peaks = 0;
            else
                # Correct the indices to account for frequency oversampling.
                peak_couples(:,1) = ...
                    (peak_couples(:,1) - 1) / oversample_factor + 1;
            endif

            # Loop through them, trying to remove one.
            removed_peak = 0;
            peak_number = 1;
            for peak_number = 1:num_peaks
                # Build analysis kernel. We make it half the section length
                # because the window function on the kernel and the window
                # function on the analysis section (constructed below) must
                # overlap 100%. If the window functions are the same size, we get
                # no time resoluton because the window functions would coincide
                # exactly.

                # We still need to subtract 1 to get the proper frequency because
                # the indexing of all these vectors starts at 1, but the spectrum
                # starts at 0 Hz.
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

                # Run section through filter.
                filter_output = fftconv(windowed_section, filter_kernel);

                # Find highest absolute peak in the response
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
#                for filt_oput_samp = 2:filter_output_length-1
#                    xm1 = abs(filter_output(filt_oput_samp-1));
#                    xp0 = abs(filter_output(filt_oput_samp));
#                    xp1 = abs(filter_output(filt_oput_samp+1));
#                    if (xm1 >= xp0 || xp0 < xp1)
#                        continue;
#                    endif
#                    c0 = filter_output(filt_oput_samp);
#                    c1 = 0.5 * (filter_output(filt_oput_samp+1) - 
#                                filter_output(filt_oput_samp-1));
#                    c2 = 0.5 * (filter_output(filt_oput_samp+1) + 
#                                filter_output(filt_oput_samp-1)) - c0;
#                    x_at_response_peak = -c1 / (2 * c2);
#                    response_peak = (c2 * x_at_response_peak + c1) * ...
#                                    x_at_response_peak + c0;
#                    if (abs(response_peak) > abs(max_response))
#                        max_response = response_peak;
#                        x_at_max_response = filt_oput_samp + x_at_response_peak;
#                        int_x_at_max_response = filt_oput_samp;
#                        frac_x_at_max_response = x_at_response_peak;
#                    endif
#                endfor
                
                # If the max response is in such a place that the corresponding
                # wavelet doesn't overlap with the windowed section 100%, this
                # peak is no good and we need to keep looking.
                salient_filt_oput_begin = kernel_length;
                salient_filt_oput_end = ...
                    salient_filt_oput_begin + window_length - kernel_length;
                if (x_at_max_response < salient_filt_oput_begin || ...
                    x_at_max_response > salient_filt_oput_end)
                   peak_number = peak_number + 1;
                   continue;
                endif
                
                # This peak is good. We will remove it.
                
                # Translate ordinates so they're based on the windowed
                # section rather than the filter output.
                int_x_max_resp_trans = int_x_at_max_response - kernel_length + 1;
                x_at_max_resp_trans = x_at_max_response - kernel_length + 1;

                # Generate wavelet based on magnitude peak.
                
                # The wavelet has the same shape as the filter kernel. We seek
                # such a scale for it that its square-sum matches the maximum
                # response value we have found. The kernel was constructed with
                # a norm of one, so we multiply the kernel's amplitude by the
                # response to get the norm of the wavelet where we want it.
                wavelet_amplitude = filter_kernel_amplitude * max_response;
                wavelet_length = kernel_length;
                
                # There's no point in letting the start magnitude of the wavelet
                # be 0.
                wavelet_t_first = -frac_x_at_max_response - wavelet_length/2 + ...
                                  ceil(frac_x_at_max_response);
                wavelet_t_last = wavelet_t_first + wavelet_length - 1;
                wavelet_window_function = ...
                    0.5 * (1 + cos([wavelet_t_first:wavelet_t_last]' * ...
                                   2 * pi / wavelet_length));
                wavelet = wavelet_amplitude * ...
                          cos([wavelet_t_first:wavelet_t_last]' * ...
                              peak_freq_rads_per_samp) .* wavelet_window_function;

                # Subtract wavelet from section
                removal_start = int_x_max_resp_trans + ...
                                ceil(frac_x_at_max_response);
                removal_end = removal_start + wavelet_length - 1;
                windowed_section(removal_start:removal_end) = ...
                    windowed_section(removal_start:removal_end) - wavelet;

                # Add wavelet to result list
                new_grain = ...
                    [x_at_max_resp_trans+section_begin-1 ...
                     peak_freq_rads_per_samp wavelet_amplitude wavelet_length];
                if (grains == 0)
                    grains = new_grain;
                else
                    grains(end+1,:) = new_grain;
                endif

                # We have removed the peak. Now see about the next one in the
                # list.
                peak_number = peak_number + 1;
                
            end
            
            # We've found all the partials to be found at this scale, so recurse.
            # This will complete the decomposition of this windowed section.
            if (window_length >= 64)
                new_grains = 0;
                [new_grains, new_residual] = ...
                    MPvR(windowed_section, window_length/2, oversample_exp, ...
                         noise_floor, depth+1);

                if (new_grains(1) ~= 0)
                    # The positions of the grains we got back are in the
                    # coordinate system used by the recursed call to MPvR.
                    # They need to be transformed to our coordinate system.
                    # The frequencies are in radians per sample, which is a
                    # global coordinate system. So is the amplitude scale.
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
            
            # The decomposition of this windowed section is complete, so subtract
            # it from the input and move on.
            corrected_input(section_begin:section_end) = ...
                corrected_input(section_begin:section_end) - ...
                section .* window_function;
                
            break;
            
        endwhile

        # We have exited from the loop that decomposes the windowed section. We
        # step to the next window. Remember that window overlap is 50%.
        
    endfor

endfunction

###############################################################################################

function plot_grains(grains, handle)
    
	grains(:,3) = log(abs(grains(:,3)));
    scf(handle);
    clf();
    plot_matrix = zeros(500,260);
    ranges = [255 499 31 255];

    for i = 1:size(grains,2);
        mins(i) = min(grains(:,i));
        maxs(i) = max(grains(:,i));
        ords(:,i) = ...
            (grains(:,i) - mins(i)) / (maxs(i) - mins(i)) * ranges(i) + 1;
    endfor
    for i = 1:size(grains,1)
        plot_matrix(501-round(ords(i,2)), [round(ords(i,1)):round(ords(i,1)+grains(i,4))-1]) = plot_matrix(501-round(ords(i,2)), [round(ords(i,1)):round(ords(i,1)+grains(i,4))-1]) + ords(i,3);
    endfor
    
    #Matplot(plot_matrix);
    #f = gcf();
    #f.color_map = bonecolormap(32);
    
endfunction

###############################################################################################

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
    spectrum_of_upsample = ...
        [spectrum(1:size(spectrum,1)/2-1) ; ...
         zeros(size(spectrum,1) * (upsample_factor-1), 1) ; ...
         spectrum(size(spectrum,1)/2:size(spectrum,1))];

    upsampled_padded_in = ifft(spectrum_of_upsample);
    out = upsampled_padded_in(1 : in_length*upsample_factor);

endfunction

###############################################################################################

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

###############################################################################################

function [c3,c2,c1,c0] = cubic_interp(xm1,xp0,xp1,xp2)
endfunction

###############################################################################################

function array_of_peaks = find_peaks(spectrum, minimum_magnitude, ...
                                     min_lobe_width, max_lobe_width)

    array_of_peaks = [0 0];
    num_peaks = 0;
    
    # Compute magnitudes
    spectrum = spectrum(:);
    mag_spectrum = abs(spectrum);

    log_minimum_magnitude = log(minimum_magnitude + 1d-300);
    mag_spectrum_length = size(mag_spectrum, 1);
    
    # Use 0 Hz to Nyquist, but also extra on the ends so peaks can
    # be found there if they exist. Nyquist is at mag_spectrum_length/2+1.
    spectrum_for_peak_search = ...
        log(1d-100 + [mag_spectrum(mag_spectrum_length+1 : mag_spectrum_length) ; ...
                      mag_spectrum(1 : mag_spectrum_length/2 + 1)]);
                                       
    spec_for_pk_search_size = size(spectrum_for_peak_search,1);
    max_peak = -1d300;
    x_at_max_peak = 0;
    peak_index = 0;
    trough_index = 1;
    found_a_lobe = 0;
#    # Take the first difference.
#    register_a = spectrum_for_peak_search(2:end) - 
#                 spectrum_for_peak_search(1:end-1);
#    register_b = sign(register_a);
#    # We can use this info to get zero crossings of the first difference.
#    register_a = register_b(2:end) - register_b(1:end-1);
#    # This will find the location of peaks.
#    register_b = find(register_a, -2);
#    # This will find the location of troughs.
#    register_c = find(register_a, 2);
#    # Distance between troughs is the lobe size.
#    register_d = register_c(2:end) - register_c(1:end-1);
#    
#
    for i = 2 : spec_for_pk_search_size-1
        # Check to see if this point is a peak, just in case later we find
        # this peak has the right width.
        if (spectrum_for_peak_search(i) > spectrum_for_peak_search(i-1) && ...
            spectrum_for_peak_search(i) >= spectrum_for_peak_search(i+1))
            peak_index = i;
        endif
        if (spectrum_for_peak_search(i) <= spectrum_for_peak_search(i-1) && ...
            spectrum_for_peak_search(i) < spectrum_for_peak_search(i+1))
            dist_from_last_trough = i - trough_index;
            trough_index = i;
            if (peak_index ~= 0)
                found_a_lobe = 1;
            endif
        endif
        if (found_a_lobe ~= 0)
            dist_from_peak = trough_index - peak_index;
            # Try to guess how wide the lobe would be if it weren't smushed
            # up against other things.
            #estimated_lobe_width = ...
            #    2 * max([dist_from_peak dist_from_last_trough-dist_from_peak]);
            # This symmetry thing is about confirming good separation of the
            # partial. max_symmetry is the reciprocal of min_symmetry because
            # we don't want to worry about the order of the arguments in the
            # division.
            min_symmetry = 0.5;   # 0.5 is an ad hoc value
            max_symmetry = 1.0 / min_symmetry;
            symmetry = dist_from_peak / (dist_from_last_trough-dist_from_peak);
            if (dist_from_last_trough < min_lobe_width || ...
                dist_from_last_trough > max_lobe_width || ...
                symmetry < min_symmetry | symmetry > max_symmetry)
               continue;
            endif
            
            # Figure out the magnitude of that peak and save it if it's
            # the maximum so far.
            c0 = spectrum_for_peak_search(peak_index);
            c1 = 0.5 * (spectrum_for_peak_search(peak_index+1) - ...
                        spectrum_for_peak_search(peak_index-1));
            c2 = 0.5 * (spectrum_for_peak_search(peak_index+1) + ...
                        spectrum_for_peak_search(peak_index-1)) - c0;
            x_at_peak = -c1 / (2 * c2);
            peak = (c2 * x_at_peak + c1) * x_at_peak + c0;
            array_of_peaks(num_peaks+1,:) = ...
                [x_at_peak+peak_index peak];
            num_peaks = size(array_of_peaks, 1);
            found_a_lobe = 0;
            peak_index = 0;
        endif
    endfor

    # Sort the peak values column, then permute the peak indices column
    # accordingly.
    [array_of_peaks(:,2), perm] = sort(array_of_peaks(:,2));
    array_of_peaks(:,1) = array_of_peaks(perm,1);

endfunction

###############################################################################################
