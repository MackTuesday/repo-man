function MPvR(in, window_length, plot_function, a, b)

in = in(:);
in_length = size(in, 1);
residual_length = ceil(in_length / (window_length/2)) * (window_length/2);
residual = [in ; zeros(residual_length-in_length, 1)];

loop_number = 0;
hop_size = window_length / 2;

    // Grab window's worth
    section_begin = 1 + loop_number * hop_size;
    section_end = section_begin + window_length - 1;
    section = residual(section_begin:section_end);

    // Apply window
    window_function = 0.5 * (1 + cos([1:window_length/2]' * ...
                                     2 * %pi / window_length));
    window_function = [ones(window_length/2,1) ; window_function];
    windowed_section = section .* window_function;

for qqq = 0:9999 do

scf(0);
clf();
plot(plot_function);
pause;
disp(qqq, 'loop number');

    // Compute energy of windowed section
    windowed_section_energy = sum(windowed_section .* windowed_section);

    // FFT
    mag_spectrum = abs(fft(windowed_section));
scf(99);
clf();
plot(log(mag_spectrum));
    // Find highest peak with parabolic interpolation, rejecting peaks that
    // aren't the right width. The main lobe of the analysis window is 4 bins
    // wide, but we allow peaks of width 3 because that's enough separation
    // for our purposes.
    // We will achieve this by taking the second difference and looking for
    // 2 or 3 negative values in a row.
    
    // Use 0 Hz to Nyquist, but also a couple extra on the ends so peaks can
    // be found there if they exist. Nyquist is at mag_spectrum_length/2+1.
    mag_spectrum_length = size(mag_spectrum, 1);
    spectrum_for_peak_search = [mag_spectrum(mag_spectrum_length-1) ; ...
                                mag_spectrum(mag_spectrum_length) ;
                                mag_spectrum(1 : mag_spectrum_length/2+3)];
    spec_for_pk_search_size = size(spectrum_for_peak_search,1);
    second_diffs = 0.25 * ...
                   spectrum_for_peak_search(1 : spec_for_pk_search_size-2);
    second_diffs = second_diffs - ...
                   0.5 * ...
                   spectrum_for_peak_search(2 : spec_for_pk_search_size-1);
    second_diffs = second_diffs +...
                   0.25 *...
                   spectrum_for_peak_search(3 : spec_for_pk_search_size);
    run_of_negs = 0;
    max_peak = 0;
    x_at_max_peak = 0;
    peak_index = 0;
    for i = 1 : mag_spectrum_length/2+1 do
        // Check to see if this point is a peak, just in case later we find
        // this peak has the right width.
        spec_for_peak_search_idx = i + 2;
        if (spectrum_for_peak_search(spec_for_peak_search_idx) > ...
            spectrum_for_peak_search(spec_for_peak_search_idx-1) & ...
            spectrum_for_peak_search(spec_for_peak_search_idx) >= ...
            spectrum_for_peak_search(spec_for_peak_search_idx+1)) then
            peak_index = spec_for_peak_search_idx;
        end
        second_diffs_index = i + 1;
        if second_diffs(second_diffs_index) < 0 then
            run_of_negs = run_of_negs + 1;
        else
            if run_of_negs == 2 | run_of_negs == 3 | run_of_negs == 4 then
                if peak_index == 0 then
                    run_of_negs = 0;
                    continue;
                end
                // Figure out the magnitude of that peak and save it if it's
                // the maximum so far.
                c0 = spectrum_for_peak_search(peak_index);
                c1 = 0.5 * (spectrum_for_peak_search(peak_index+1) - ...
                            spectrum_for_peak_search(peak_index-1));
                c2 = 0.5 * (spectrum_for_peak_search(peak_index+1) + ...
                            spectrum_for_peak_search(peak_index-1)) - c0;
                x_at_peak = -c1 / (2 * c2);
                peak = (c2 * x_at_peak + c1) * x_at_peak + c0;
                if peak > max_peak then
                    max_peak = peak;
                    x_at_max_peak = peak_index + x_at_peak;
                end
            else
                peak_index = 0;
            end
            run_of_negs = 0;
        end
    end
    
    // A final check for a run is required.
    if run_of_negs == 2 | run_of_negs == 3 | run_of_negs == 4 then
        // Figure out the magnitude of that peak and save it if it's
        // the maximum so far.
        c0 = spectrum_for_peak_search(peak_index);
        c1 = 0.5 * (spectrum_for_peak_search(peak_index+1) - ...
                    spectrum_for_peak_search(peak_index-1));
        c2 = 0.5 * (spectrum_for_peak_search(peak_index+1) + ...
                    spectrum_for_peak_search(peak_index-1)) - c0;
        x_at_peak = -c1 / (2 * c2);
        peak = (c2 * x_at_peak + c1) * x_at_peak + c0;
        if peak > max_peak then
            max_peak = peak;
            x_at_max_peak = peak_index + x_at_peak;
        end
    end
    
    // If no peak is found, analyze section of half size.
    if x_at_max_peak == 0 then
disp(window_length/2, '*** recursing ***');
scf(99);
clf();
plot(log(mag_spectrum));
//error('Debug quit');
        if window_length >= 16 then
            result_section = 0;
            
            temp_section = windowed_section(1:window_length/2);
            result_section = MPvR(temp_section, window_length/2, ...
                                  plot_function, a, a+window_length/2-1);
            temp_section = temp_section - result_section;
            windowed_section(1:window_length/2) = temp_section;
            temp_section = plot_function(a:a+window_length/2-1);
            temp_section = temp_section - result_section;
            plot_function(a:a+window_length/2-1) = temp_section;

            temp_section = ...
                windowed_section(window_length/4+1:window_length*3/4);
            result_section = MPvR(temp_section, window_length/2, ...
                                  plot_function, a+window_length/4, ...
                                  a+window_length*3/4-1);
            temp_section = temp_section - result_section;
            windowed_section(window_length/4+1:window_length*3/4) = ...
                temp_section;
            temp_section = ...
                plot_function(a+window_length/4:a+window_length*3/4-1);
            temp_section = temp_section - result_section;
            plot_function(a+window_length/4:a+window_length*3/4-1) = ...
                temp_section;

            temp_section = ...
                windowed_section(window_length/2+1:window_length);
            result_section = MPvR(temp_section, window_length/2, ...
                                  plot_function, a+window_length/2, ...
                                  a+window_length-1);
            temp_section = temp_section - result_section;
            windowed_section(window_length/2+1:window_length) = ...
                temp_section;
        end
    
        return windowed_section;
    end
    
    // Build analysis kernel. We make it half the section length because the
    // window function on the kernel and the window function on the analysis
    // section (constructed below) must overlap 100%. If the window functions
    // are the same size, we get no time resoluton because the window
    // functions would coincide exactly.
    // max_peak was an index into spectrum_for_peak_search.
    // Correct it to point into mag_spectrum.
    x_at_max_peak = x_at_max_peak - 2;
    // We still need to subtract 1 to get the proper frequency because the
    // indexing of all these vectors starts at 1, but the spectrum starts at
    // 0 Hz.
disp(x_at_max_peak - 1, 'peak_freq');
    peak_freq_rads_per_samp = 2 * %pi * (x_at_max_peak - 1) / ...
                              mag_spectrum_length;
    kernel_length = window_length / 2;
    kernel_window_function = ...
        0.5 * (1 + cos([-kernel_length/2+1:kernel_length/2]' * ...
                       2 * %pi / kernel_length));
    filter_kernel = cos([-kernel_length/2+1:kernel_length/2]' * ...
                        peak_freq_rads_per_samp) .* ...
                    kernel_window_function;
    filter_kernel_amplitude = 1 / sqrt(sum(filter_kernel .* filter_kernel))
    filter_kernel = filter_kernel * filter_kernel_amplitude;

    // Grab analysis section
//    analysis_section_begin = section_begin;
//    analysis_section_end = analysis_section_begin + window_length - 1;
//    analysis_section = residual(analysis_section_begin:analysis_section_end);
    analysis_section = windowed_section;

    // Run analysis section through filter.
    filter_output = conv(analysis_section, filter_kernel);
    extended_filter_output = filter_output;
    //extended_filter_output = complex_extender(filter_output);
//scf(94);
//clf();
//plot(filter_output);
    // Complex-extend filter output
    // We care about only the portion of the filter output that corresponds
    // to 100% overlap between the filter kernel and the analysis section.
    salient_filt_oput_begin = kernel_length;
    salient_filt_oput_end = salient_filt_oput_begin + window_length - ...
                            kernel_length;
    salient_filter_output = ...
        extended_filter_output(salient_filt_oput_begin:salient_filt_oput_end);

    // Find highest peak in magnitude.
    salient_filt_oput_mags = abs(salient_filter_output);
    salient_filt_oput_len = size(salient_filt_oput_mags, 1);
    max_response = 0;
    x_at_max_response = 0;
    frac_x_at_max_response = 0;
    for i = 2:salient_filt_oput_len-1
        if (salient_filt_oput_mags(i) <= salient_filt_oput_mags(i-1) | ...
            salient_filt_oput_mags(i) <  salient_filt_oput_mags(i+1)) then
            continue;
        end
        c0 = salient_filt_oput_mags(i);
        c1 = 0.5 * (salient_filt_oput_mags(i+1) - ...
                    salient_filt_oput_mags(i-1));
        c2 = 0.5 * (salient_filt_oput_mags(i+1) + ...
                    salient_filt_oput_mags(i-1)) - c0;
        x_at_response_peak = -c1 / (2 * c2);
        response_peak = (c2 * x_at_peak + c1) * x_at_peak + c0;
        if response_peak > max_response then
            max_response = response_peak;
            x_at_max_response = i + x_at_response_peak;
            int_x_at_max_response = i;
            frac_x_at_max_response = x_at_response_peak;
        end
    end

    if max_response == 0 then
        x_at_max_response = 1;
        frac_x_at_max_response = 0;
        max_response = salient_filt_oput_mags(1);
        int_x_at_max_response = 1;
        if salient_filt_oput_mags(salient_filt_oput_len) > max_response then
            max_response = salient_filt_oput_mags(salient_filt_oput_len);
            int_x_at_max_response = salient_filt_oput_len;
        end
    end

    if salient_filter_output(int_x_at_max_response) < 0 then
        max_response = -max_response;
    end

disp(max_response, 'max response');

    // Measure the energy of the input around the peak so we know how to
    // scale the wavelet we're about to generate.
    energy_estimate_begin = int_x_at_max_response - 1;
    if energy_estimate_begin < 1 then
        energy_estimate_begin = 1;
    end
    energy_estimate_end = int_x_at_max_response - 1 + kernel_length - 1;
    analysis_sect_around_pk = ...
       analysis_section(energy_estimate_begin:energy_estimate_end);
    em1 = sqrt(sum(analysis_sect_around_pk .* analysis_sect_around_pk));
    
    energy_estimate_begin = int_x_at_max_response;
    energy_estimate_end = int_x_at_max_response - 1 + kernel_length - 1 + 1;
    analysis_sect_around_pk = ...
       analysis_section(energy_estimate_begin:energy_estimate_end);
    ep0 = sqrt(sum(analysis_sect_around_pk .* analysis_sect_around_pk));
    
    energy_estimate_begin = int_x_at_max_response + 1;
    energy_estimate_end = int_x_at_max_response - 1 + kernel_length - 1 + 2;
    if energy_estimate_end > size(analysis_section,1) then
        energy_estimate_end = size(analysis_section, 1);
    end
    analysis_sect_around_pk = ...
       analysis_section(energy_estimate_begin:energy_estimate_end);
    ep1 = sqrt(sum(analysis_sect_around_pk .* analysis_sect_around_pk));
    
    c0 = ep0;
    c1 = 0.5 * (ep1 - em1);
    c2 = 0.5 * (ep1 + em1) - c0;
    est_energy_at_the_peak = (c2 * frac_x_at_max_response + c1) * ...
                             frac_x_at_max_response + c0;

    // Generate wavelet based on magnitude peak
    // There's no point in letting the start magnitude of the wavelet be 0.
//    wavelet_amplitude = filter_kernel_amplitude * est_energy_at_the_peak * ...
//                        sign(max_response);
    wavelet_amplitude = filter_kernel_amplitude * max_response;
disp(wavelet_amplitude, 'wavelet_amplitude');
if wavelet_amplitude >= 1 then
    pause;
end
    wavelet_length = kernel_length;
    wavelet_sample_first = -frac_x_at_max_response - wavelet_length/2 + 1;
    wavelet_sample_last = wavelet_sample_first + wavelet_length - 1;
    wavelet_window_function = ...
        0.5 * (1 + cos([wavelet_sample_first:wavelet_sample_last]' * ...
                       2 * %pi / wavelet_length));
    wavelet = wavelet_amplitude * ...
              cos([wavelet_sample_first:wavelet_sample_last]' * ...
                  peak_freq_rads_per_samp) .* wavelet_window_function;

scf(95);
clf();
plot(conv(windowed_section, filter_kernel));

    // Subtract wavelet from section
    removal_start = round(x_at_max_response) + 1;
    removal_end = removal_start + wavelet_length - 1;
    windowed_section(removal_start:removal_end) = ...
        windowed_section(removal_start:removal_end) - wavelet;

    temp_function = plot_function(a+removal_start-1:a+removal_end-1);
    temp_function = temp_function - wavelet;
    plot_function(a+removal_start-1:a+removal_end-1) = temp_function;
    
padded_wavelet = [zeros(removal_start-1,1) ; wavelet ; ...
                  zeros(window_length-removal_end,1)];
scf(96);
clf();
plot(conv(padded_wavelet, filter_kernel));
q=5;
    // Compute new section energy

end

// Remove section from residual

// Shift location of next window

endfunction
