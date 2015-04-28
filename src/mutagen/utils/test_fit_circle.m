function [xs,ys,r2s] = test_fit_circle(S)

    S = S(:);
    slen = size(S,1);
    
    for i = 1:slen-2
        [xs(i),ys(i),r2s(i)] = fit_circle([real(S(i));imag(S(i))], ...
                                          [real(S(i+1));imag(S(i+1))], ...
                                          [real(S(i+2));imag(S(i+2))]);                               
    endfor

endfunction