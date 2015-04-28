function sgramplot2(y)

  cmaprows = 64;
  mycmap = zeros(cmaprows*cmaprows, 3);
  mycmap(end-cmaprows+1:end,:) = hsv(cmaprows);
  for i = 1:cmaprows-1
    istart = 1 + cmaprows * (i-1);
    iend = istart + cmaprows - 1;
    mycmap(istart:iend, :) = mycmap(end-cmaprows+1:end, :) * i / cmaprows;
  end
  
  myimgmag = log(abs(y) + 1e-300) / log(10);
  % The 3 corresponds to 10^-3 = -60 dB
  myimgmag = floor((myimgmag + 3) / 3 * cmaprows) + 1;
  myimgmag = max(myimgmag, 1);
  myimgmag = min(myimgmag, cmaprows);
  
  myimgarg = floor(mod(arg(y) + 2*pi, 2*pi) * cmaprows / (2*pi));
  myimgarg = max(myimgarg, 1);
  myimgarg = min(myimgarg, cmaprows);

  myimg = myimgarg + (myimgmag-1) * cmaprows;
  colormap(mycmap);  
  image(myimg);
#  image([1:cmaprows] + [0:cmaprows:cmaprows*(cmaprows-1)]');
  
end