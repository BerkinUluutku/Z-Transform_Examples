function [fitresult,gof] = springpot_fit(freq_axis_rad_cut,Q_cut_abs, lnr)
[xData, yData] = prepareCurveData( freq_axis_rad_cut, Q_cut_abs );

% Set up fittype and options.
fit_str = ...
    sprintf('abs(C * ((1-exp(-%.4f)*exp(-1i*x))^b))',lnr);
ft = fittype( 'abs(C * ((1-exp(-0.075)*exp(-1i*x))^b))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.MaxFunEvals = 10000;
opts.MaxIter = 10000;
opts.StartPoint = [360000000 0.022];
opts.TolFun = 1e-20;
opts.TolX = 1e-20;
opts.Upper = [Inf 1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end

