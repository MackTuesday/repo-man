function y = exp_window(N, r1, r2)

  domain = [0:N-1]' / N;
  w1 = exp(domain * -r1);
  w1 = w1 - w1(end);
  w2 = 1.0 - exp(domain * -r2);
  w2 = w2 ./ w2(end);
  y = w1 .* w2;
  y /= max(y);

endfunction