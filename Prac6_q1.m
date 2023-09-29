clear
clf

hold on
[DayNum, TideHeight] = ReadWATideData('WYWYN01_04_2019.txt', 2019);
plot(DayNum, TideHeight, 'b-')
PeriodDays = GetTidalConstituentPeriods(5);

AMat = BuildTidalLSQCoefftMat(DayNum, PeriodDays);

stddev = 0.16;
YVec = TideHeight;
N = length(TideHeight);
WMat = speye(N) / stddev^2;
M1 = (AMat' * WMat * AMat);
M2 = (AMat' * WMat * YVec);
ThetaVec = M1 \ M2;
ModelVec = AMat * ThetaVec;

plot(DayNum, ModelVec, 'r-')


ResVec = ModelVec - TideHeight;
plot(DayNum, ResVec, 'g-')
legend('Measurements', 'Model', 'Residuals', fontsize=16)

chi2 = ResVec' * WMat * ResVec;

unitVar = chi2/(length(TideHeight) - length(ThetaVec));
hold off
disp(chi2)
disp(unitVar)

CovarMat = unitVar * inv(AMat' * WMat * AMat);
varMat = diag(CovarMat);
stdDevMat = sqrt(varMat);

disp(std(ResVec))




