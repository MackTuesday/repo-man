function y = matreduce(x, rpows, cpows)

  rfactor = 2 ^ rpows;
  cfactor = 2 ^ cpows;
  y = zeros(ceil(rows(x)/rfactor), ceil(columns(x)/cfactor));
  
  for r = 1:rfactor
    for c = 1:cfactor
      sumrows = ceil((rows(x)-r+1)/rfactor);
      sumcols = ceil((columns(x)-c+1)/cfactor);
      y(1:sumrows, 1:sumcols) = y(1:sumrows, 1:sumcols) + x(r:rfactor:end, c:cfactor:end);
    endfor
  endfor
  
  y = y / (rfactor * cfactor);

endfunction