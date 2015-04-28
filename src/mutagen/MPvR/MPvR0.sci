function MPvR(in, window_length)

in = in(:);
in_length = size(in, 1);
residual = [zeros(window_length/2-1,1) ; in ; zeros(window_length/2-1,1)];

loop_number = 0;
hop_size = window_length / 2;

// Grab window's worth
// Note that because of the zero padding and the progression of whittling
// these sections to nothing, the first section will be half zeros.
section_begin = 1 + loop_number * hop_size;
section_end = section_begin + window_length - 1;
section = residual(section_begin, section_end);

// Apply window
window_function = 
    [ones(window_length/2-1, 1)' ; ...
     0.5*(1+cos([0:window_length/2-1] * 2 * %pi * window_length))];
section = section .* window_function;

// Quit if window energy is too small
threshold = 0.001;
section_energy = sum(section .* section);
if section_energy < threshold then
    return;
end

while section_energy >= threshold do

    // FFT
    section_spectrum = fft(section);

    // Find highest peak with parabolic interpolation, but skip the last
    // half of the spectrum because it's redundant
    max_peak = 0;
    x_at_max_peak = 0;
    for i = 2:section_length-1
        if (section(i) <= section(i-1) | section(i) < section(i+1)) then
            peaks(i) = 0;
            continue;
        end
        c0 = section(i);
        c1 = 0.5 * (section(i+1) - section(i-1));
        c2 = 0.5 * (section(i+1) + section(i-1));
        x_at_peak = -c1 / (2 * c2);
        peak = ((c2 * x_at_peak) + c1) * x_at_peak) + c0;
        if peak > max_peak then
            max_peak = peak;
            x_at_max_peak = i + x_at_peak;
        end
    end

    // Build analysis kernel
    // Only window_length-1 support samples are needed because the peak of
    // the kernel is in the exact center. Also, remember that convolution
    // performs dot products with the reverse of the filter kernel. It doesn't
    // matter here because cosine.
    filter_kernel = cos([-window_length/2+1:window_length/2-1] * ...
                        2 * %pi * (x_at_max_peak - 1) / window_length) .* ...
                    window_function;

    // Grab analysis section
    // The first half of the section is zeros. The peak response can be
    // anywhere in the section interval. So put a half window interval
    // through the filter.
    analysis_section_begin = section_begin + window_length/2;
    analysis_section_end = analysis_section_begin + window_length/2 - 1;
    analysis_section = residual(analysis_section_begin, analysis_section_end);

    // Run analysis section through filter. Take only the middle half-window's
    // worth of the filter output because that corresponds to the interval of
    // the analysis section that we really care about.
    filter_output = conv(analysis_section, filter_kernel);
    filter_output_middle_start = window_length;
    filter_output_middle_end = filter_output_middle_start + window_length/2);
    filter_output = filter_output(window_length:window_length);
    filter_output = filter_output(:);
    return;

    // Complex-extend filter output
    extended_filter_output = complex_extender(filter_output);

    // Find highest peak in magnitude
    filter_output_mags = abs(extended_filter_output);
    filter_output_length = size(filter_output, 1);
    max_peak = 0;
    x_at_max_peak = 0;
    frac_x_at_max_peak = 0;
    peaks(:) = 0;
    for i = 2:filter_output_length-1
        if (filter_output_mags(i) <= filter_output_mags(i-1) | ...
            filter_output_mags(i) <  filter_output_mags(i+1)) then
            peaks(i) = 0;
            continue;
        end
        c0 = filter_output_mags(i);
        c1 = 0.5 * (filter_output_mags(i+1) - filter_output_mags(i-1));
        c2 = 0.5 * (filter_output_mags(i+1) + filter_output_mags(i-1));
        x_at_peak = -c1 / (2 * c2);
        peak = ((c2 * x_at_peak) + c1) * x_at_peak) + c0;
        if peak > max_peak then
            max_peak = peak;
            x_at_max_peak = i + x_at_peak;
            frac_x_at_max_peak = x_at_peak;
        end
    end

    // Generate wavelet based on magnitude peak
    // There's no point in letting the start magnitude of the wavelet be 0.
    // If the end magnitude is 0 it's because frac_x_at_max_peak was exactly 0.
//    wavelet_sample_first = -frac_x_at_max_peak - window_length/2 + 1;
//    wavelet_sample_last = wavelet_sample_first + window_length - 1;
//    wavelet = max_peak * cos([wavelet_sample_first:wavelet_sample_last] * ...
//                             2 * %pi * (x_at_max_peak-1) / window_length) .* ...
//              window_function;

    // Subtract wavelet from input to get residual
//    wavelet_residual_start = 
//    residual(wavelet_residual_start:wavelet_residual_end) = ...
//      residual(wavelet_residual_start:wavelet_residual_end) - wavelet;
//    section(wavelet_section_start:wavelet_section_end) = ...
//      section(wavelet_section_start:wavelet_section_end) - ...
//      wavelet(wavelet_section_overlap_start:wavelet_section_overlap_end);

    //section_energy = sum(section .* section);
//    section_energy = 0;

end

// Shift location of next window
endfunction
