function sgramplot(y)

  cmapsize = 256;
  mycmap = jet(cmapsize);
  myimg = log(y + 1e-300) / log(10);
  % The 3 corresponds to 10^-3 = -60 dB
  myimg = floor((myimg + 3) / 3 * cmapsize) + 1;  
  myimg = max(myimg, 1);  
  myimg = min(myimg, cmapsize);  
  colormap(mycmap);  
  image(myimg);  
  
end