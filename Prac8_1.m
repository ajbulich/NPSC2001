function Prac8_1()
%Function to demo minimisation of a function of many parameters.
%
%The problem is as follows:
%Consider N identical masses of mass M, connected together by identical
%springs with unstretched length L0 and spring constant Ks.  The two end
%masses are connected to the top of vertical posts of  heights H1 and H2
%which are separated horizontally by distance XDist.
%
%Write a program that can find the equilibrium positions of the masses when given
%specific values of N, M, L0, Ks, H1 and H2 by minimising the potential
%energy of the system.
%
%Alec Duncan, Curtin University, April 2020

%First save all the parameters that define the problem in a structure so we
%can easily pass it through to the cost function

FigNum = 11;
PlotIntermedSolns = false;

Info.H1 = 4.0;  %Height of first post, (m)
Info.H2 = 2.0;  %Height of second post (m)
Info.H3 = 3.0;
Info.H4 = 3.5;
Info.XDist = 3.0;  %Horizontal separation between posts 1 and 2 (m)
Info.YDist = 4.0; %Horizontal separation between posts 2 and 3 (m)
Info.XDist2 = 3.5; %Horizontal separation between posts 3 and 4 (m)
Info.YDist2 = 5.0; %Horizontal separation between posts 4 and 1 (m)
Info.NMass = 4^2; %Number of masses
Info.Mass = 2.0/Info.NMass;  %Mass of each (kg)

%Calculate the unstretched length as a proportion of the straight line
%distance between the end points (note there are N+1 springs)
SLDist = sqrt(Info.XDist^2 + (Info.H2 - Info.H1)^2);
SLDist2 = sqrt(Info.YDist^2 + (Info.H3 - Info.H2)^2);
SLDist3 = sqrt(Info.XDist2^2 + (Info.H4 - Info.H3)^2);
SLDist4 = sqrt(Info.YDist2^2 + (Info.H1 - Info.H4)^2);
Info.L0 = 0.9 * SLDist/(sqrt(Info.NMass+1));
Info.L0_2 = 0.9 * SLDist2/(sqrt(Info.NMass+1));  %added sqrt because the total masses
Info.L0_3 = 0.9 * SLDist3/(sqrt(Info.NMass+1));  %is the masses in one direction squared
Info.L0_4 = 0.9 * SLDist4/(sqrt(Info.NMass+1));
Info.L0x = mean([Info.L0, Info.L0_3]);
Info.L0y = mean([Info.L0_2, Info.L0_4]);
Info.Ks = 3*Info.NMass;  %Spring constant (N/m)
Info.NMass2 = sqrt(Info.NMass);

Info.g = 9.8;  %Acceleration due to gravity (m/s^2)


%Coordinate system is X horizontal, Y Horizontal, Z positive up

%Set up the vector of initial parameters, this will be all the X values
%followed by all the Z values in one big column vector.  Started with the masses on a straight line
%between the end points.
dX = Info.XDist/sqrt(Info.NMass + 1);
dY = Info.YDist/sqrt(Info.NMass + 1);
X0 = (1:Info.NMass2).' * dX;
Y0 = (1:Info.NMass2).' * dY;
Z0 = (1:Info.NMass).' * (Info.H2 - Info.H1)/(Info.NMass + 1) + Info.H1;

[xmat, ymat] = meshgrid(X0, Y0);
BigXVec = xmat(:);
BigYVec = ymat(:);
BigZVec = Z0;

Par0 = [BigXVec; BigYVec; BigZVec]; % 16 X's then 16 Y's then 16 Z's
disp(Par0)



XParMin = BigXVec - Info.XDist/10;
XParMax = BigXVec + Info.XDist/10;
YParMin = BigYVec - Info.YDist/10;
YParMax = BigYVec + Info.YDist/10;
ZParMin = zeros(size(BigZVec));
ZParMax = BigZVec + 0.5;

ParMin = [XParMin; YParMin; ZParMin];
ParMax = [XParMax; YParMax; ZParMax];


NPar = length(Par0);
XPosts = [0; Info.XDist; Info.XDist2; 0];
YPosts = [0; 0; Info.YDist; Info.YDist2];
ZPosts = [Info.H1; Info.H2; Info.H3; Info.H4];
postVec = [XPosts; YPosts; ZPosts];
hold on
plot3(Par0(1:NPar/3), Par0(1+NPar/3:2*NPar/3), Par0(1+2*NPar/3:end), 'ro')
plot3(XPosts, YPosts, ZPosts, 'bo')
hold off
%Plot this up to make sure we got it right

Opt.NewtonDamping = 1;
Opt.NReset = 500;  %Conjugate gradient is reset every this many iterations
Opt.AbsTol = 1e-7;
Opt.RelTol = 1e-7;
Opt.MaxIt = 20000;
%Handle to function that calculates the gradient of the cost function
%(Jacobian).  Empty to use numerical derivatives
%Opt.GradFnHndl = @lCatenaryGradFn;  
Opt.GradFnHndl = [];

%Method = {'Newton'};
disp(' ');
disp('Method, Seconds taken, Number of iterations, Number of cost fn calls, Number of gradient fn calls, Final energy')
%Method = {'Gradient descent', 'Conjugate gradient', 'Newton', 'BFGS'};
%Method = {'Conjugate gradient', 'BFGS'};
Method = {'Conjugate gradient'};

for IMeth = 1:length(Method)
    figure(FigNum);
    FigNum = FigNum + 1;
    clf;
    subplot(2,1,1);
    hold on;
    %Plot the end points
    plot3([0, Info.XDist], [0, Info.YDist], [Info.H1, Info.H2],[Info.XDist, Info.XDist2], [Info.YDist, Info.YDist2], [Info.H2, Info.H3],[Info.XDist2, Info.XDist], [Info.YDist2, Info.YDist], [Info.H3, Info.H4], 'ro', 'MarkerFaceColor', 'r');
    lPlotPar(Par0, Info);
    lCatenaryCostFn('Reset');
    if ~isempty(Opt.GradFnHndl)
        Opt.GradFnHndl('Reset');
    end
    tic;
    switch Method{IMeth}
        case 'Reduced gradient'
            Opt.UseHessian = true;
            [ParVec, FnVal, ParamTrace, FnTrace] = LibNDMinimiseRG(Info, [], Par0, @lRGDSCatCostFn, ParMin, ParMax, false, Opt);
        case 'Direction set'
            [ParVec, FnVal, ParamTrace, FnTrace] = LibNDMinimise(Info, [], Par0, @lRGDSCatCostFn, ParMin, ParMax, false, Opt);
            
        case 'Simulated annealing'
            Opt.TFactor = 0.99;
            [ParVec, FnVal, ParamTrace, TempTrace, FnTrace, AllParamTrace, AllFnTrace] = ...
                LibSimulatedAnneal(Info, [], Par0, @lRGDSCatCostFn, XParMin, XParMax, false, Opt);
            %ParamTrace = AllParamTrace;
        otherwise
            
            [ParVec, FnVal, ParamTrace, FnTrace, RetStatus] =  UnconstrainedMin(@lCatenaryCostFn, Par0, Method{IMeth}, Info, Opt);
    end
    ElapsedTime = toc;
    if PlotIntermedSolns
        for ITrace = 1:length(FnTrace)
            lPlotPar(ParamTrace(:, ITrace), Info);
        end
    end
    
    if ~isempty(ParVec)
        plot3(ParVec(1:NPar/3), ParVec(NPar/3+1:2*NPar/3), ParVec(2*NPar/3+1:end), 'ko', 'MarkerFaceColor', 'k')
    end
    
    NIt = length(FnTrace);

    title([Method{IMeth} ', ' num2str(Info.NMass) ' masses']);
    xlabel('X (m)');
    ylabel('Z (m)');
    axis([0 Info.XDist 0 1.1*max(Info.H1, Info.H2)]);
    grid on;
    box on;
    
    subplot(2,1,2);
    plot(FnTrace);
    xlabel('Iteration');
    ylabel('Potential energy (J)');
    %axis([0 inf 20 30]);
    grid on;
    box on;
    drawnow;
    title([num2str(NIt) ' iterations, final energy = ' num2str(FnVal) ' J']);
    CostFnCount = lCatenaryCostFn('GetCount');
    if ~isempty(Opt.GradFnHndl)
        GradFnCount = Opt.GradFnHndl('GetCount');
    else
        GradFnCount = 0;
    end

    disp([Method{IMeth} ', ' num2str(ElapsedTime) ', ' num2str(NIt) ', ' num2str(CostFnCount) ', ' num2str(GradFnCount) ', ' num2str(FnVal)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lPlotPar(ParVec, Info)
NPar = length(ParVec);
XPlt = [0; ParVec(1:NPar/3); Info.XDist];
YPlt = [0; ParVec(NPar/3+1:2*NPar/3); Info.YDist];
ZPlt = [Info.H1; ParVec(2*NPar/3+1:end); Info.H2];
plot3(XPlt, YPlt, ZPlt, '.-');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PE = lCatenaryCostFn(ParVec, Info)
%Stuff so I can count the times this is called
persistent CallCount;
if ischar(ParVec)
    switch ParVec
        case 'Reset'
            CallCount = 0;
            PE = [];
        case 'GetCount'
            PE = CallCount;
    end
else
    %Actual cost function calcs
    CallCount = CallCount + 1;
    
    NPar = length(ParVec);
    XVec = ParVec(1:NPar/3);
    YVec = ParVec(NPar/3+1:2*NPar/3);
    ZVec = ParVec(2*NPar/3+1:end);
    
    %Total gravitational potential energy
    PEGravity = sum(Info.Mass * Info.g * ZVec);
    
    dX = [diff(XVec(1:Info.NMass2))];
    dX2 = [diff(XVec(Info.NMass2*(Info.NMass2-1)+1:end));];
    dX3 = [diff(XVec(Info.NMass2+1:0.5*(Info.NMass2)^2));];
    dY = [diff(YVec(1:Info.NMass2));];
    dY2 = [diff(YVec(Info.NMass2*(Info.NMass2-1)+1:end));];
    dY3 = [diff(YVec(Info.NMass2+1:0.5*(Info.NMass2)^2));];
    dZ = [diff(ZVec(1:Info.NMass2));];
    dZ2 = [diff(ZVec(Info.NMass2*(Info.NMass2-1)+1:end));];
    dZ3 = [diff(ZVec(Info.NMass2+1:0.5*(Info.NMass2)^2));];
    
    Dist1 = sqrt(dX.^2 + dZ.^2 + dY.^2);
    Dist2 = sqrt(dX2.^2 + dY2.^2 + dZ2.^2);
    Dist3 = sqrt(dX3.^2 + dY3.^2 + dZ3.^2);
    PESprings_x = sum(0.5*Info.Ks*(Dist1 - Info.L0x).^2) + (Info.NMass2-2)*sum(0.5*Info.Ks*(Dist3 - Info.L0x).^2) + sum(0.5*Info.Ks*(Dist2 - Info.L0x).^2);
    PESprings_y = sum(0.5*Info.Ks*(Dist1 - Info.L0y).^2) + (Info.NMass2-2)*sum(0.5*Info.Ks*(Dist3 - Info.L0y).^2) + sum(0.5*Info.Ks*(Dist2 - Info.L0y).^2); 
    PESprings_posts = 0.5*Info.Ks*(norm([0-XVec(1),0-YVec(1),Info.H1-YVec(1)]))^2 + 0.5*Info.Ks*(norm([Info.XDist-XVec(Info.NMass2),0-YVec(Info.NMass2),Info.H2-ZVec(Info.NMass2)]))^2 + 0.5*Info.Ks*(norm([Info.XDist2-XVec(2*Info.NMass2+1),Info.YDist-YVec(2*Info.NMass2+1),Info.H3-ZVec(2*Info.NMass2+1)]))^2 + 0.5*Info.Ks*(norm([0-XVec(Info.NMass2*(Info.NMass2-1)+1),Info.YDist2-YVec(Info.NMass2*(Info.NMass2-1)+1),Info.H4-ZVec(Info.NMass2*(Info.NMass2-1)+1)]))^2;
    PESprings = PESprings_y + PESprings_x + PESprings_posts;
    PE = PEGravity + PESprings;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FnVal = lRGDSCatCostFn(ExtraInfo, Z, Params)
%This provides a slightly different interface to the cost funciton for
%compatibility with some pre-existing optimisation routines
FnVal = lCatenaryCostFn(Params, ExtraInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Grad = lCatenaryGradFn(ParVec, Info)
%Analytic calculation of the gradient of the cost function
%To avoid approximations in the numerical calc and speed things up

%Stuff so I can count the times this is called
persistent CallCount;
if ischar(ParVec)
    switch ParVec
        case 'Reset'
            CallCount = 0;
            Grad = [];
        case 'GetCount'
            Grad = CallCount;
    end
else
    %Actual gradient function calcs
    CallCount = CallCount + 1;
    
    NPar = length(ParVec);
    XVec = ParVec(1:NPar/2);
    ZVec = ParVec(NPar/2+1:end);
    
    XAll = [0; XVec; Info.XDist];
    ZAll = [Info.H1; ZVec; Info.H2];
    
    dX = diff(XAll);
    dZ = diff(ZAll);
    
    DistAll = sqrt(dX.^2 + dZ.^2);
    Fact = 1 - (Info.L0 ./ DistAll);
    
    dUdX = - Info.Ks * (Fact(2:end) .* dX(2:end) - Fact(1:end-1) .* dX(1:end-1));
    
    dUdZ = - Info.Ks * (Fact(2:end) .* dZ(2:end) - Fact(1:end-1) .* dZ(1:end-1)) + ...
        Info.Mass * Info.g;
    
    Grad = [dUdX; dUdZ];
end
