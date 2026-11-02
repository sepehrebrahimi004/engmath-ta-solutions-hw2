%MAIN  Engineering Mathematics - TA solutions, MATLAB homework 2: the Fourier
%      transform (symbolic and FFT), convolution, scaling, image phase and
%      musical-note identification.

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

%% Q1  Duality, modulation and a cosine filter (symbolic, frequency f in Hz)
syms t f real
FT = @(x) simplify(subs(fourier(x, t, sym('w')), sym('w'), 2*sym(pi)*f));   % X(f) = int x e^{-j2pi f t}
h  = sinc(t);  g = rectangularPulse(t);
H  = FT(h);    G = FT(g);
fprintf('Q1  F{sinc(t)} = %s,   F{rect(t)} = %s\n', char(H), char(G));
fig = figure('Visible', 'off', 'Position', [100 100 1100 600]);
subplot(2, 3, 1); fplot(h, [-6 6]); grid on; title('h(t) = sinc(t)');
subplot(2, 3, 4); fplot(H, [-2 2]); grid on; title('H(f) = rect(f)'); xlabel('f (Hz)'); ylim([-0.1 1.1]);
subplot(2, 3, 2); fplot(g, [-2 2]); grid on; title('g(t) = rect(t)'); ylim([-0.1 1.1]);
subplot(2, 3, 5); fplot(G, [-6 6]); grid on; title('G(f) = sinc(f)'); xlabel('f (Hz)');
% modulation: cos(2 pi t) h(t) moves H to +-1 Hz at half height
H1 = FT(cos(2*pi*t)*h);
subplot(2, 3, 3); fplot(H, [-3 3], 'k'); hold on; fplot(H1, [-3 3], 'r', 'LineWidth', 1.2); grid on;
ylim([-0.1 1.1]); title('F\{cos(2\pit) sinc t\}: shifted to \pm1 Hz'); legend('H', 'H_1');
% multiplying by cos(5 omega) = cos(10 pi f) gives two echoes: [h(t+5) + h(t-5)]/2
h2 = (subs(h, t, t + 5) + subs(h, t, t - 5))/2;
H2 = FT(h2);
fprintf('    check: F{[h(t+5)+h(t-5)]/2} - cos(10 pi f) H(f) = %s\n', char(simplify(rewrite(H2 - cos(10*pi*f)*H, 'exp'))));
subplot(2, 3, 6); fplot(h, [-10 10], 'k'); hold on; fplot(h2, [-10 10], 'r', 'LineWidth', 1.2); grid on;
title('cos(5\omega)H(\omega) \rightarrow [h(t+5) + h(t-5)]/2'); legend('h', 'h_2');
save_figure(fig, 'q1_duality_modulation');

%% Q2  rect * rect = triangle, so F{triangle} = sinc^2  (FFT with proper scaling)
dt = 1e-3;  tt = -4:dt:4;
r  = double(abs(tt) < 0.5);
tri = conv(r, r, 'same')*dt;                               % the dt makes the discrete sum an integral
Nf = 2^nextpow2(numel(tt)*8);                              % zero padding for a smooth spectrum
fx = (-Nf/2:Nf/2-1)/(Nf*dt);
R  = fftshift(abs(fft(r, Nf)))*dt;   T = fftshift(abs(fft(tri, Nf)))*dt;
band = abs(fx) <= 5;
fprintf('Q2  triangle peak %.3f (expected 1);  max |FFT - sinc| = %.1e,  max |FFT - sinc^2| = %.1e\n', ...
        max(tri), max(abs(R(band) - abs(sinc(fx(band))))), max(abs(T(band) - sinc(fx(band)).^2)));
fig = figure('Visible', 'off', 'Position', [100 100 1000 550]);
subplot(2, 2, [1 2]); plot(tt, r, tt, tri, 'LineWidth', 1.3); grid on; xlim([-2 2]);
legend('rect(t)', 'rect * rect = \Lambda(t)'); xlabel('t');
subplot(2, 2, 3); plot(fx, R, 'b', fx, abs(sinc(fx)), 'r--', 'LineWidth', 1.1); grid on; xlim([-5 5]);
title('|FFT\{rect\}| vs |sinc(f)|'); xlabel('f (Hz)');
subplot(2, 2, 4); plot(fx, T, 'b', fx, sinc(fx).^2, 'r--', 'LineWidth', 1.1); grid on; xlim([-5 5]);
title('FFT\{\Lambda\} vs sinc^2(f)'); xlabel('f (Hz)');
save_figure(fig, 'q2_convolution_theorem');

%% Q3  Time scaling: F{rect(a t)} = (1/|a|) sinc(f/a)
a = [0.2 0.5 1 2 4];
fig = figure('Visible', 'off', 'Position', [100 100 1000 420]);
for k = 1:numel(a)
    Fa = FT(rectangularPulse(a(k)*t));
    subplot(1, 2, 1); fplot(rectangularPulse(a(k)*t), [-4 4], 'LineWidth', 1.2); hold on;
    subplot(1, 2, 2); fplot(Fa, [-6 6], 'LineWidth', 1.2); hold on;
    fprintf('Q3  a = %-4g  F(0) = %-5g (1/a = %g)   F(f = a) = %g (first zero at f = a)\n', ...
            a(k), double(limit(Fa, f, 0)), 1/a(k), double(subs(Fa, f, a(k))));
end
subplot(1, 2, 1); grid on; title('rect(a t)'); ylim([-0.1 1.1]);
subplot(1, 2, 2); grid on; title('(1/a) sinc(f/a)'); xlabel('f (Hz)');
legend(arrayfun(@(x) sprintf('a = %g', x), a, 'UniformOutput', false));
save_figure(fig, 'q3_time_scaling');

%% Q4  Swapping FFT magnitude and phase of two colour images (per channel)
A = im2double(imresize(imread('peppers.png'), [384 512]));          % sample images shipped with MATLAB
B = im2double(imresize(imread('ngc6543a.jpg'), [384 512]));
FA = fft2(A);  FB = fft2(B);                                          % fft2 acts on each colour channel
AB = real(ifft2(abs(FA) .* exp(1j*angle(FB))));                       % |A|, phase of B
BA = real(ifft2(abs(FB) .* exp(1j*angle(FA))));                       % |B|, phase of A
fig = figure('Visible', 'off', 'Position', [100 100 1000 620]);
subplot(2, 2, 1); imshow(A); title('A');
subplot(2, 2, 2); imshow(B); title('B');
subplot(2, 2, 3); imshow(rescale(AB)); title('|A| with the phase of B');
subplot(2, 2, 4); imshow(rescale(BA)); title('|B| with the phase of A');
save_figure(fig, 'q4_phase_swap');
fprintf('Q4  corr(|B| + phase A, A) = %.2f,  corr(|B| + phase A, B) = %.2f\n', ...
        corr(BA(:), A(:)), corr(BA(:), B(:)));

%% Q5  Identify seven notes from their spectra
% The course recordings are replaced by synthesised notes: a decaying fundamental plus
% harmonics, at two sample rates and in a shuffled order, like the original set.
names = {'C4 (Do)', 'D4 (Re)', 'E4 (Mi)', 'F4 (Fa)', 'G4 (Sol)', 'A4 (La)', 'B4 (Si)'};
freqs = 440 * 2.^(([-9 -7 -5 -4 -2 0 2])/12);
order = [4 5 1 6 2 7 3];                                    % hidden sequence F G C A D B E
rates = [44100 44100 44100 44100 44100 48000 44100];
harm  = [1 0.6 0.35 0.2 0.1];
rng(7);
fig = figure('Visible', 'off', 'Position', [100 100 1100 750]);
fprintf('Q5  %-6s %-12s %-10s %s\n', 'note', 'peak (Hz)', 'detected', 'true');
for k = 1:7
    fs = rates(k);  tn = 0:1/fs:1.3;  f0 = freqs(order(k)) * (1 + 0.002*randn);   % slightly out of tune
    x = sum(harm.' .* sin(2*pi*(1:5).'*f0*tn), 1) .* exp(-2.5*tn) + 0.01*randn(size(tn));
    Nf = 2^18;  X = abs(fft(x, Nf))/numel(x);  fa = (0:Nf-1)*fs/Nf;
    keep = fa < 1500;  [~, i] = max(X(keep));  fpk = fa(i);
    [~, j] = min(abs(12*log2(fpk./freqs)));                 % nearest equal-tempered pitch
    fprintf('    %-6d %-12.1f %-10s %s\n', k, fpk, names{j}, names{order(k)});
    subplot(4, 2, k); plot(fa(keep), X(keep), 'LineWidth', 1); grid on; xlim([0 1500]);
    title(sprintf('note %d: peak %.1f Hz \\rightarrow %s', k, fpk, names{j}));
end
subplot(4, 2, 7); xlabel('f (Hz)'); subplot(4, 2, 8); axis off;
text(0, 0.5, {'Each note shows its fundamental and', 'harmonics at 2f_0, 3f_0, ...;', ...
     'the pitch is the nearest semitone to the strongest peak.'});
save_figure(fig, 'q5_note_spectra');
