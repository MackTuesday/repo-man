function Y = complex_extender(X)
    X = X(:);
    X_length = size(X, 1);
    padded_X_length = 0;
    for i = 0:31
        p = 2^i;
        if p >= X_length then
            padded_X_length = p;
            break;
        end
    end
    
    if padded_X_length == 0 then
        error('You suck');
        return [0];
    end
    
    padded_X = [X ; zeros(padded_X_length - X_length, 1)];
    X_spectrum = fft(padded_X);
    extended_X_spectrum = X_spectrum;
    // H(w) = F(w)*s(w)
    // s = i for w>0; -i for w<0, 0 for w=0
    // i*-i(a + bi) + (a + bi) =  a+bi + a+bi = 2a+2bi ==> a+bi = F(w)
    // i* i(a - bi) + (a - bi) = -a+bi + a-bi = 0
    extended_X_spectrum(padded_X_length/2+1:padded_X_length) = 0
    extended_X = fft(extended_X_spectrum, 1);
    Y = extended_X(1:X_length);
endfunction
