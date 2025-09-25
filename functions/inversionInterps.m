eps_real = real(eps_model(:));
eps_imag = imag(eps_model(:));
% Hessian Interp (Occurd Within Inversion Debug)
x = eps_real(:);
y = eps_imag(:);
z = J_weightNorm(:);
[eps_unique, ~, ic] = unique([x y], 'rows');
z_avg = accumarray(ic, z, [], @mean);
JweightInterp = scatteredInterpolant(eps_unique(:,1), eps_unique(:,2), z_avg, 'natural', 'nearest');

% Regularization Interp (Occured after parameter optimization)
lambdaJ_flat = lambdaJ(:);
lambdaT_flat = lambdaT(:);
[eps_unique, ~, ic] = unique([eps_real eps_imag], 'rows');
lambdaJ_avg = accumarray(ic, lambdaJ_flat, [], @mean);
lambdaT_avg = accumarray(ic, lambdaT_flat, [], @mean);
lambdaJInterp = scatteredInterpolant(eps_unique(:,1), eps_unique(:,2), lambdaJ_avg, 'natural', 'nearest');
lambdaTInterp = scatteredInterpolant(eps_unique(:,1), eps_unique(:,2), lambdaT_avg, 'natural', 'nearest');