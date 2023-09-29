clear


T = 0.1;      % sinusoidal period
N = 10;     % number of pulses
TSig = T*N; % total signal length
TZeroes = 22;
TMax = TSig + TZeroes;% total length including zeros
if mod(TMax,2) ~= 0
    TMax = TMax + 1;
end
f0 = 1/T;   % sinusoid frequency
A = 2;      % sinusoid amplitude

dT = 0.01;  % Sample interval (s)
Fs = 1/dT;  % Sample rate (Samples/s)

TVec = 0:dT:TMax-dT;  %Build time vector
NPt = length(TVec);  %Number of points in total
NSig = TSig/TMax * NPt;

Sig = A*sin(2*pi*f0*TVec);
plot(TVec, Sig);
for ii=1:NPt
    if ii>NSig
        Sig(ii) = 0;
    end
end

plot(TVec, Sig)

Spec = fft(Sig);

dF = Fs/NPt;  %Frequency interval
FPos = (0:NPt/2)*dF;  %Positive frequencies

FNeg = -fliplr(FPos(2:end-1)); 
FVec = [FPos FNeg]; % non flipped version
size(FVec)

FPlt = [FNeg FPos];  % flipped version
size(FPlt)
SpecPlt = [Spec(NPt/2+2:end) Spec(1:NPt/2+1)];
SpecPlt = 1/max(SpecPlt) * SpecPlt;


fDom = -24.999:0.001:25;
for f2=1:1:50000
    f = (f2-25000)/1000;
    E(f2) = (0.5*T*N*sin(N*pi*T*(f-f0)))/(pi*T*N*(f-f0)) - (0.5*T*N*sin(N*pi*T*(f+f0)))/(pi*N*T*(f+f0));
end
E = 1/max(E) * E;

figure(2);
lPlotSpectrum(FPlt, SpecPlt, 'Re-ordered',fDom, E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotSpectrum(FVec, Spec, TitleStr, fDom, E)

clf;
subplot(2,1,1);
hold on;
plot(FVec, real(Spec));
plot(FVec, imag(Spec));
plot(FVec, abs(Spec));
plot(fDom, E);
legend('Real', 'Imaginary', 'Maginitude', 'Analytic');
ylabel('Spectrum');
grid on;
box on;
title(TitleStr);


subplot(2,1,2);
plot(FVec, rad2deg(angle(Spec)));
ylabel('Phase (deg)');
xlabel('Frequency (Hz)');
grid on;
box on;

end

