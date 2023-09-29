clear
[Y, FS] = audioread('Wintergatan - Marble Machine (music instrument using 2000 marbles).mp3', [48000*20 48000*50]);
figure(1)
plot(Y(:, 1))
title("Waveform")
xlabel('Samples')
NFFT = 4000;
NOverlap = 0.5*NFFT; 

[S, F, T] = spectrogram(Y(:,1), hanning(NFFT), NOverlap, NFFT, FS);

C = S .* conj(S);
C = 20*log10(C);
figure(2)
pcolor(T, F, C)
shading interp;
colorbar;
Hnd = colorbar;
Hnd.Label.String = 'dB';
set(gca, 'YScale', 'log')
clim([-50 50]);
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')



[Pxx, F2] = pwelch(Y(:, 1), hanning(NFFT), NOverlap, NFFT, FS);

figure(3)
plot(F2, Pxx, 'b-')
ylabel('Power Density')
xlabel('Frequency (Hz)')
set(gca, 'XScale', 'log')
