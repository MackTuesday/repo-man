function y = spectrogram(x, w, overlapnumber)

  x = x(:);
  w = w(:);
  
  xlen = size(x, 1);
  winlength = size(w, 1);
  stride = winlength / overlapnumber;
  numstrides = (floor(xlen/stride) + overlapnumber - 1);
  fullxlen = numstrides * stride;
  fullx = [x ; zeros(fullxlen-xlen, 1)];

  y = zeros(winlength, numstrides);
  
  for winstart = 1:stride:(fullxlen-winlength)  
    n = (winstart - 1) / stride + 1;  
    winend = winstart + winlength - 1;  
    fftin = w .* fullx(winstart:winend);
    y(:, n) = abs(fft(fftin));  
  endfor  
  
  y = y / winlength; 

endfunction