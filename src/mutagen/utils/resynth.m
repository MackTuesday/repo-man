function out = resynth(sg, overlap, fadein, fadeout)

    numrows = size(sg, 1);
    numcols = size(sg, 2);
    
    stride = numrows / overlap;
    totallength = stride * (numcols + overlap - 1);
    out = zeros(totallength, 1);
    fadeinlength = floor(numrows * fadein);
    fadeinpart = [0:fadeinlength-1]' / fadeinlength;
    fadeoutlength = floor(numrows * fadeout);
    fadeoutpart = [fadeoutlength:-1:1]' / fadeoutlength;
    w = [fadeinpart; ones(numrows-fadeinlength-fadeoutlength,1); fadeoutpart];

    for i = 1:numcols
        addstart = (i-1) * stride + 1;
        addend = addstart + numrows - 1;
        out(addstart:addend) = out(addstart:addend) + w .* real(ifft(sg(:,i)));
    endfor

endfunction