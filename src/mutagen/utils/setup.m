more off
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\wavesmono.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\introx.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether A.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether A L.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether A R.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether A C.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether B.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether B L.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether B R.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\Tether B C.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\YCtL.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\YCtL L.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\YCtL R.wav')';
X = wavread('c:\\users\\owner\\dropbox\\mutagen\\YCtL C.wav')';
X = wavread('c:\\users\\Brent\\dropbox\\mutagen\\Tether A C.wav')';
w = asym_window(0, 512, 1, 3, 0, 0);
sg = stftbrent(X, w, 4);
softsg = hack(sg, 0, 7, 0.25, 0);
softout = resynth(softsg, 4, 0.05, 0.15);
wavwrite(softout/max(abs(softout)), 44100, 16, "c:\\users\\owner\\desktop\\softout.wav");
multisgram(X, w, 4, 6, 0, 4);
image(floor(rand(10,10)*100))
colormap("hot")
colormap("jet")
colormap("winter")
colormap("ocean")
map = autumn(256);

X = wavread('c:\\users\\owner\\dropbox\\mutagen\\wavesmono.wav')';
S = sin([0:8191] * 2 * pi / 64)';
WS = cos([0:8191] * 2 * pi / 2048)' .* 0.5 .* (1 - cos([0:8191] * 2 * pi / 8192))';
more off
[gr,res] = MPvR(S, 4096);
[gr,res] = MPvR(S, 4096, 0, 1.0e-10, 0);
PS = sin([0:8191] * 2 * pi / 64)' + sin([0:8191] * 2 * pi / 49)' + sin([0:8191] * 2 * pi / 81)';
PSfrac = sin([0:8191] * 2 * pi / 64.32)' + sin([0:8191] * 2 * pi / 49.57)' + sin([0:8191] * 2 * pi / 81.19)';
PSphase = sin([9:8200] * 2 * pi / 64)' + sin([3:8194] * 2 * pi / 49)' + sin([0:8191] * 2 * pi / 81)';
S49 = sin([0:8191] * 2 * pi / 49)';
[gr,res] = MPvR(PS, 4096, 0, 1.0e-10, 0);
Sf = sin([0:8191] * 2 * pi / 82.84)';
[gr,res] = MPvR(Sf, 4096, 0, 1.0e-10, 0);
B = cos([-255:255]' * 2 * pi / 256) .* 0.5 .* (1 + cos([-255:255]' * 2 * pi / 512));
B = [zeros(129,1) ; B ; zeros(129,1)];
[gr,res] = MPvR(B, 8192, 0, 1.0e-10, 0);
w = 0.5 * (1.0 - cos([0:8191] * 2 * pi / 8192));
dw_dt = -imag(ifft(fft(w) .* [8191:-1:0] * 2 * pi / 8192));
qqq = fr_spectrum(PSfrac, w, dw_dt, 3);


S25 = sin([0:8191] * 2 * pi / 25)';
S31 = sin([0:8191] * 2 * pi / 31)';
S49 = sin([0:8191] * 2 * pi / 49)';
S64 = sin([0:8191] * 2 * pi / 64)';
S81 = sin([0:8191] * 2 * pi / 81)';
S121 = sin([0:8191] * 2 * pi / 121)';
S169 = sin([0:8191] * 2 * pi / 169)';
S289 = sin([0:8191] * 2 * pi / 289)';
S361 = sin([0:8191] * 2 * pi / 361)';
Na = 2 * rand(8192,1) - 1;
Nb = 2 * rand(8192,1) - 1;
Nc = 2 * rand(8192,1) - 1;
PS = S31 + S81 + S289;
predictors = [S25 S31 S49 S64 S81 S121 S169 S289 S361 Na Nb Nc];


Pd = dlmread("c:\\users\\brent\\dropbox\\mutagen\\prostate.data", "\t", 1, 1);
Pd(:,10) = [];
y = Pd(:,9);
Pd(:,9) = [];
F = fopen("c:\\users\\brent\\dropbox\\mutagen\\prostate.data");
titles = fgetl(F);
titles = regexprep(titles, '\s+', ' ');
col = strsplit(titles, ' ');
inputLabel = cellstr(col);
inputLabel(11) = [];
outputLabel = col(10);
inputLabel(10) = [];
inputLabel(1) = [];
fclose(F);
[c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10] = textreadx('c:\\users\\brent\\dropbox\\mutagen\\prostate.data', '%d %d %d %d %d %d %d %d %d %d %s', 'headerlines', 1);
train = char(c10);
idx = (train == 'T');
y = y(idx);
X = Pd(idx,:);