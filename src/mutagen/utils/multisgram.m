function multisgram(x, w, overlapnumber, numscales, freduce, treduce)

  x = x(:);
  w = w(:);
  ypartsizes = zeros(numscales, 1);
  y = [];
  
  for scalenumber = 1:numscales
    winsubsample = 2 ^ (scalenumber - 1);
    win = w(1:winsubsample:end);
    sg = stftbrent(x, win, overlapnumber);
    ypart = matreduce(sg(2:end/2+1,:), freduce, treduce);
%figure(scalenumber);
%sg = sg / max(max(sg));
%sgramplot(sg);
    ypartrows(scalenumber) = size(ypart, 1);
    ypartcols(scalenumber) = size(ypart, 2);
    ypartsize = ypartrows(scalenumber) * ypartcols(scalenumber);
    y = [y ; reshape(ypart, ypartsize, 1)];
  endfor
  
  mmax = abs(max(max(y)));
  y = y / mmax;
  ystart = 1;
  for scalenumber = 1:numscales
    ypartsize = ypartrows(scalenumber) * ypartcols(scalenumber);
    yend = ystart + ypartsize - 1;
    ypart = reshape(y(ystart:yend), ypartrows(scalenumber), ypartcols(scalenumber));
    figure(scalenumber);
    sgramplot2(ypart(end:-1:1, :));
    ystart = yend + 1;
  endfor

endfunction