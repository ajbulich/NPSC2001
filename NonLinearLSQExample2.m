function NonLinearLSQExample2()
%This example simulates measured data using the
%equation y = ax + b*exp(c*x), adds noise to it,
%and then uses the simple iterative method to carry out a nonlinear
%least squares fit of the coefficients a, b and c
%
%Alec Duncan,  Curtin University, 15/4/2019

NRepeat = 1;  %Number of times to repeat experiment

Exact_a = 1.0;
Exact_b = 0.1;
Exact_c = 0.5;

%Simulation parameters
SigmaNoise = 1.5;  %Noise standard deviation
XVec = -5:1:10;

%Initial guess values of parameters in same order:
a0 = 0.9;
b0 = 0.12;
c0 = 0.4;

% a0 = 0.5;
% b0 = 0.2;
% c0 = 0.1;

ConvLim = 0.001;  %Stop when the norm of the change in the parameter vector is less than this
MaxIterations = 1000;  %If it hasn't converged in this many iterations, abort

Theta_Exact = [Exact_a; Exact_b; Exact_c];

%Simulate some data
YExact = lMyFun(XVec, Theta_Exact);
NMeas = length(YExact);

for IRepeat = 1:NRepeat
    
    YMeas = YExact + randn(1, NMeas) * SigmaNoise;
    
    %Plot the simulated exact and "measured" data points
    
    if NRepeat == 1
        figure(1);
        clf;
        hold on;
        plot(XVec, YExact, 'r+');
        plot(XVec, YMeas, 'bx');
    end
    
    %Set up the weight matrix.  W doesn't depend on Theta0 so this can be done
    %before the loop
    WMat = speye(NMeas) / (SigmaNoise^2);
    
    %Do the iterative least squares fit
    Theta0 = [a0; b0; c0];
    
    NPar = length(Theta0);
    
    if IRepeat == 1
        FitParMat = zeros(NRepeat, NPar);  %Matrix to store all fitted parameters
    end
    
    Converged = false;
    Aborted = false;
    Count = 0;
    
    while ~Converged && ~Aborted
        %Set up
        YfVec =  YMeas.' - lMyFun(XVec, Theta0).';
        AfMat = lMyJacobian(XVec, Theta0);
        
        %Do the least squares solution
        DeltaVec = (AfMat.' * WMat * AfMat)\(AfMat.' * WMat * YfVec);
        
        %Update the assumed parameter vector
        Theta0 = Theta0 + DeltaVec;
        
        %Optional extras:
        %Count the iterations
        Count = Count + 1;
        %Calculate the unit variance on each iteration
        ResidVec = (lMyFun(XVec, Theta0) - YMeas).';
        UnitVar = ResidVec.' * WMat * ResidVec / (NMeas - length(Theta0));
        NormDelta = norm(DeltaVec);
        disp(['Iteration: ' num2str(Count) ', norm(Delta): ' num2str(NormDelta) ...
            ', Unit variance: ' num2str(UnitVar)]);
        
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


        
        
        YFitted = lMyFun(XVec, Theta0);
        if NRepeat == 1
            plot(XVec, YFitted, 'go');
            legend('True', 'Measured', 'Fitted');
            xlabel('X')
            ylabel('Y');
            grid on;
            
            figure(2);
            clf;
            plot(XVec, ResidVec, 'x');
            grid on;
            xlabel('X')
            ylabel('Residuals')
        end
    end
    if NRepeat > 1
        CalcCovMat = cov(FitParMat);
        disp('Covariance matrix from repeated runs:')
        lDisplayCovMat(CalcCovMat);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YVals = lMyFun(XVals, ThetaVec)
%This eveluates the function we are trying to fit
%Change the function to  fit something different
YVals = ThetaVec(1)*XVals + ThetaVec(2) * exp(ThetaVec(3) * XVals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AfMat = lMyJacobian(XVals, ThetaVec)
%Function to evaluate the Jacobian matrix
%Preferably using analytic derivatives
%This also needs changing if you change the function to be fitted.

NPt = length(XVals);
Col1 = XVals.';
Col2 = exp(ThetaVec(3)*XVals.');
Col3 = XVals.' * ThetaVec(2) .* exp(ThetaVec(3)*XVals.');
AfMat = [Col1 Col2 Col3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
