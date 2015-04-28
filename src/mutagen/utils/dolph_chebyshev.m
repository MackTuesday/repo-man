function y = dolph_chebyshev(length, floor_level_db)

    floor_level = 10^(floor_level_db/20);
    x0 = cosh(acosh(1/floor_level)/length)
    W = zeros(length,1);
    x = x0 * cos(pi * [-length/2+1:length/2]' / length);
    for i = 1:length
        if (abs(x(i)) > 1)
            W(i) = cosh(length * acosh(x(i))) * floor_level;
        else
            W(i) = cos(length * acos(x(i))) * floor_level;
        endif
    endfor
    W = [W(1:length/2+1) ; zeros(size(W,1)-length,1) ; W(length/2+2:length)];
    y = abs(ifft(real(W)));
    y = y / max(y);

endfunction