function out = error_sequence(input, basis, coeffs)

    num_coeffs = size(coeffs, 2);
    out = zeros(num_coeffs, 1);
    
    for i = 1:num_coeffs
        erri = input - basis * coeffs(i,:)';
        out(i) = sqrt(sum(erri .* erri));
    endfor
    
endfunction