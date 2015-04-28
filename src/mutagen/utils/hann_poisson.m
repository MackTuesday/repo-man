function out = hann_poisson(length)

    % Force length to be even
    length = ceil(length/2) * 2;
    
    % Pre-allocate output array
    out = zeros(length, 1);
    
    slope_end = length/2;
    func_domain = [-slope_end:slope_end-1]';
    out = 0.5 * (1 + cos(func_domain * 2 * pi / (length+1))) .* ...
          exp(-32 * abs(func_domain) / (length+1));
    
endfunction
