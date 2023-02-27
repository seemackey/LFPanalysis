function [MyCoeffs,jit_tol] = taufxn_v3(x,y,tauguess,rangeguess,makecurve)


% inputs are the data and guesses for the first two parameters

xData=x';
yData=y;

% settings for the fit
ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [rangeguess tauguess 0.0001];
opts.Upper = [1 1 1];
opts.StartPoint = [0.01 0.01 0.01];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% clf;
% A = fitresult;
% h = plot( fitresult, xData, yData );

% gets coefficients 
MyCoeffs = coeffvalues(fitresult);

% make a curve so we can do activities with it
if makecurve == 1
    curvex = min(xData):0.0001:max(xData);
    curvey = MyCoeffs(1)*exp(-curvex/MyCoeffs(2))+MyCoeffs(3);
    
    % find the first place the curve has insignificant y vals
    if min(curvey) > 0.15
        jit_tol(1) = max(xData);
        jit_tol(2) = NaN;
    elseif max(yData) < 0.15
        jit_tol(1) = NaN;
        jit_tol(2) = 1;
    else
        below  = curvey<0.15;
        belowidx = find(below,1);
        jit_tol(1) = curvex(belowidx);
        jit_tol(2) = NaN;
    end
    %
    
end

end