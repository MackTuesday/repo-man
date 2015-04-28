
function corrections = fr_spectrum(in, w, oversample_exp)
    
    in = in(:);

    in_length = size(in, 1);
    oversample_factor = 2^oversample_exp;
    
    dw_dt = -imag(ifft(fft(w) .* [in_length-1:-1:0]' * 2 * pi / in_length));
    
    windowed_input = in(:) .* w(:);
    diff_windowed_input = in(:) .* dw_dt(:);
    
    padded_WI = [windowed_input ; zeros((oversample_factor-1)*in_length, 1)];
    padded_DWI = [diff_windowed_input ; zeros((oversample_factor-1)*in_length, 1)];
    
    spectrum = fft(padded_WI);
    spectrum_diff = fft(padded_DWI);

    
    spectrum_conj = conj(spectrum);
    unnormed_corrections = -imag(spectrum_diff .* spectrum_conj);
    corrections = unnormed_corrections ./ (spectrum .* spectrum_conj);

endfunction
