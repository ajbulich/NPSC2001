p1 = [-1; 2];
p2 = [2; 3];
p3 = [-0.5; -2];
p4 = [5; 0];
Points = {p1; p2; p3; p4};
XPoints = [-1 2 -0.5 5];
YPoints = [2 3 -2 0];
x_exact = -5;
y_exact = -5;
x_guess = 75000;
y_guess = 166000;
ThetaExact = [x_exact; y_exact];
Theta0 = [x_guess; y_guess];
NPoints = 4; 
NRepeat = 100;
NPar = 2;
NMeas = 4;
ConvLim = 0.001;  %Stop when the norm of the change in the parameter vector is less than this
MaxIterations = 1000;

d1 = norm(p1 - ThetaExact);
d2 = norm(p2 - ThetaExact);
d3 = norm(p3 - ThetaExact);
d4 = norm(p4 - ThetaExact);

SigmaNoise = 0.2;
i = 0;

for IRepeat = 1:NRepeat
    d1Noise = d1 + randn(1, 1)*SigmaNoise;
    d2Noise = d2 + randn(1, 1)*SigmaNoise;
    d3Noise = d3 + randn(1, 1)*SigmaNoise;
    d4Noise = d4 + randn(1, 1)*SigmaNoise;
    
    YMeas = [d1Noise, d2Noise, d3Noise, d4Noise]';
    
    WMat = eye(NPoints)/(SigmaNoise)^2;
    
    FitParMat = zeros(NRepeat, NPar);
    
    Converged = false;
    Aborted = false;
    Count = 0;
    
    while ~Converged && ~Aborted
        YfVec =  YMeas - lMyFun2(XPoints, YPoints, Theta0);
        AfMat = lMyJacobian2(XPoints, YPoints, Theta0);
    
        %Do the least squares solution
        DeltaVec = (AfMat.' * WMat * AfMat)\(AfMat.' * WMat * YfVec);
            
        %Update the assumed parameter vector
        Theta0 = Theta0 + DeltaVec;
            
        %Optional extras:
        %Count the iterations
        Count = Count + 1;
        %Calculate the unit variance on each iteration
        ResidVec = (lMyFun2(XPoints, YPoints, Theta0) - YMeas);
        UnitVar = ResidVec.' * WMat * ResidVec / (NMeas - length(Theta0));
        NormDelta = norm(DeltaVec);
        
    
        %Check for convergence and divergence
        if NormDelta < ConvLim
            Converged = true;
        else
            if Count > MaxIterations  %Need to make sure it does actually converge!
                Aborted = true;
            end
        end
    end
    if Aborted
        disp(['Failed to converge in ' num2str(Count) ' iterations']);
    else
        disp(['Number of iterations: ' num2str(Count)]);
        disp('Fitted parameters:');
        disp(Theta0);
        FitParMat(IRepeat, :) = Theta0.';
            
        ParamCovMat = inv(AfMat.' * WMat * AfMat);
        disp('Estimated paramater covariance matrix: ')
        lDisplayCovMat(ParamCovMat);
      
        YFitted = lMyFun2(XPoints, YPoints, Theta0);
       
    end
    if NRepeat > 1
        CalcCovMat = cov(FitParMat);
        disp('Covariance matrix from repeated runs:')
        lDisplayCovMat(CalcCovMat);
        disp(Theta0)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normVec2 = lMyFun2(XPoints, YPoints, Theta0)
%This eveluates the function we are trying to fit
%Change the function to  fit something different
normVec2 = [];
for j=1:4
    normVec2(j) = norm(Theta0 - [XPoints(j); YPoints(j)]);
end
normVec2 = normVec2';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AfMat = lMyJacobian2(XPoints, YPoints, Theta0)
%Function to evaluate the Jacobian matrix
%Preferably using analytic derivatives
%This also needs changing if you change the function to be fitted.
normVec = [];
for i=1:4
    normVec(i) = norm(Theta0 - [XPoints(i); YPoints(i)]);
end
normVec = normVec';
Col1_1 = (-1*XPoints+Theta0(1))';
Col1 = Col1_1 ./ normVec;
Col2_1 = (-1*YPoints+Theta0(2))';
Col2 = Col2_1 ./ normVec;
AfMat = [Col1 Col2];
end

function lDisplayCovMat(CovMat)
%Displays covariance matrix info in a few different ways

disp(CovMat)

ParamSDev = sqrt(diag(CovMat));
disp('Parameter standard deviations:')
disp(ParamSDev)

%Calculate the parameter correlation matrix - this has elements
%Gamma(i,j) = Sigma(i,j)/sqrt(Var(i) * Var(j))
%Could do it in a loop but its easier to use meshgrid
[SDev_i, SDev_j] = meshgrid(ParamSDev, ParamSDev);
ParamCorrelMat = CovMat ./ (SDev_i .* SDev_j);
disp('Parameter correlation matrix:')
disp(ParamCorrelMat)
end

