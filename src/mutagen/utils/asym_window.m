function y = asym_window(T, N, M1, M2, Zpre, Zpost)

  % T is for type. It isn't used yet.
  
  winlength1 = N * M1;
  winlength2 = N * M2;
  domain1 = [0:winlength1-1]' / winlength1;
  domain2 = [0:winlength2-1]' / winlength2;
  w1 = domain1;
  w2 = 1 - domain2;
%  w1 = 0.5 * (1.0 - cos(domain1 * pi));
%  w2 = 0.5 * (1.0 + cos(domain2 * pi));
  
  zeroslength1 = N * Zpre;
  zeroslength2 = N * Zpost;
  z1 = zeros(zeroslength1, 1);
  z2 = zeros(zeroslength2, 1);

  y = [z1 ; w1 ; w2 ; z2];

endfunction