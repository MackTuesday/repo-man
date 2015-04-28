function out = hack(X, type, x1, x2, x3, x4, x5, x6)

    numrows = size(X, 1);
    numcols = size(X, 2);
    
    out = zeros(numrows, numcols);
    for i = 1:numcols
        out(:, i) = scalepeaks(X(:, i), x1, x2, x3);
    endfor

endfunction