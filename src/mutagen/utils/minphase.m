function y = minphase(x)

    x = x(:);
    
    x = [x; zeros(size(x,1),1)];
    n = size(x,1);
    lm = log(abs(fft(x))+1e-300);     % Fourier of cepstrum -- symmetric because fft is conjugate
                                      %  symmetric
    w = [ones(n/2,1); -ones(n/2,1)];  % Hilbert thingy for flipping over acausal components
    curve = imag(fft(w.*ifft(lm)));   % minimized phase response -- antisymmetric because fft is
                                      %  conjugate symmetric because
                                      %  ifft(lm) is real because lm is symmetric
    y = ifft(exp(lm+i*curve));        % lm+i*curve is conjugate symmetric because curve is antisymmetric
    y = y(1:n/2);                     % take first half because we zero-padded original time domain

endfunction