function out = lorentz(length)

    % Force length to be even
    length = ceil(length/2) * 2;
    
    % Pre-allocate output array
    out = zeros(length, 1);
    
    half_length = length/2;
    func_half_preimage = [0:half_length]' / half_length;
    func_preimage = [func_half_preimage(2:half_length+1) ; ...
                     func_half_preimage(half_length:-1:1)];
    out = 1 - sqrt(1 - func_preimage.*func_preimage);
    
endfunction
