clear
a = 0.2;
d = 1;
N = 20;

Length_One = a + d;
Length_Total = Length_One * N;

dx = 0.005;  % Sample interval (m)
us = 1/dx;  % Sample rate (Samples/m)
NPt_One = Length_One / dx;
NPt = NPt_One * N;
aPt = 0.2/Length_One * NPt_One;
dPt = NPt_One - aPt;

% This loops creates one signal pulse
for i=1:NPt_One
    if i<aPt
        Sig(i) = 1;
    else
        Sig(i) = 0;
    end
end
% Now to write code to create N of them
Sig1 = Sig;
for ii = 2:N
    Sig = [Sig Sig1];
end

plot(1:NPt, Sig)
% perfect, now to do fourier transform and all that lovely stuff.



Spec = fft(Sig);

du = us/NPt;  %Frequency interval
FPos = (0:NPt/2)*du;  %Positive frequencies

FNeg = -fliplr(FPos(2:end-1)); 
FVec = [FPos FNeg]; % non flipped version
size(FVec)

FPlt = [FNeg FPos];  % flipped version
size(FPlt)
SpecPlt = [Spec(NPt/2+2:end) Spec(1:NPt/2+1)];
size(SpecPlt)

SpecPlt = 1/max(SpecPlt) * SpecPlt;
figure(2);

uDom = -24.99:0.01:25;
for u2=1:1:5000
    u = (u2-2500)/100;
    sum = 0;
    for n=1:N
        sum = sum + exp(2*pi*1i*u*(n-1)*(a+d));
    end
    num = abs(sum);
    E(u2) = a*(sin(pi*u*a))/(pi*u*a) .* num ;
end
E = 1/max(E) * E;


lPlotSpectrum(FPlt, SpecPlt, 'Re-ordered', uDom, E);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotSpectrum(FVec, Spec, TitleStr, uDom, E)

clf;
subplot(2,1,1);
hold on;
plot(FVec, real(Spec));
plot(FVec, imag(Spec));
plot(FVec, abs(Spec));
plot(uDom, E);
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

