function lwc_extended = extrapolateLWC(lwc_data, n_extra, method)
    % lwc_data: original LWC vector (1D)
    % n_extra: number of points to extrapolate
    % method: 'exp' or 'power'

    N = numel(lwc_data);
    x = (1:N)';
    y = lwc_data(:);

    % Define extrapolation model
    switch lower(method)
        case 'exp'
            modelFun = @(b,x) b(1) * exp(-b(2) * x) + b(3);
            beta0 = [max(y), 0.01, min(y)];
        case 'power'
            modelFun = @(b,x) b(1) * x.^b(2);
            beta0 = [1, -1];
    end

    % Fit
    options = optimset('Display','off');
    beta = lsqcurvefit(modelFun, beta0, x, y, [], [], options);

    % Extrapolate
    x_ext = (N+1:N+n_extra)';
    y_ext = modelFun(beta, x_ext);

    % Add synthetic noise (decaying magnitude)
    % decay_sigma = std(y(end-100:end)) * 0.75;
    % noise = decay_sigma * randn(size(x_ext));

    % Final result
    % lwc_extended = [y; max(6e-5, y_ext + noise)];  % clamp negative
    % lwc_extended = [max(6e-5, y_ext + noise)];  % clamp negative
    lwc_extended = [max(6e-5, y_ext)];  % clamp negative
end