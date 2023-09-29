clear
clf

hold on
[DayNum, TideHeight] = ReadWATideData('WYWYN01_04_2019.txt', 2019);
[DayNum2, TideHeight2] = ReadWATideData('WYWYN01_04_2020.txt', 2019);
plot(DayNum2, TideHeight2, 'b-')
PeriodDays = GetTidalConstituentPeriods(5);

AMat = BuildTidalLSQCoefftMat(DayNum, PeriodDays);
AMat2 = BuildTidalLSQCoefftMat(DayNum2, PeriodDays);

stddev = 0.16;
YVec = TideHeight;
N = length(TideHeight);
WMat = speye(N) / stddev^2;
M1 = (AMat' * WMat * AMat);
M2 = (AMat' * WMat * YVec);
ThetaVec = M1 \ M2;
ModelVec = AMat * ThetaVec;
ModelVec2 = AMat2 * ThetaVec;


plot(DayNum2, ModelVec2, 'r-')


ResVec = ModelVec2 - TideHeight2;
plot(DayNum2, ResVec, 'g-')
legend('2020 Measurements', '2019 Model', 'Residuals', fontsize=16)

disp(std(ResVec))




