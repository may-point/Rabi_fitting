data =readtable('rabi_fixed_freq1.xlsx','Sheet','sheet_1')
tau = data{:,1};       % tau in nanoseconds
A_data = data{:,2};    % amplitude

% Plot scatter
figure;
scatter(tau, A_data, 'filled'); hold on;
xlabel('\tau (ns)');
ylabel('Amplitude');
title('Rabi Oscillations with Damping');

% Define model: A = Asin(2*pi*Omega_R*tau + phi) * exp(-tau/T2) + C
damped_sine = @(b, tau) b(1)*sin(2*pi*b(2)*tau + b(3)) .* exp(-tau / b(4)) + b(5);

% Initial guesses
b0 = [1.6, 0.0083, 1.57, 350, -3.2];  % [A, Ω_R (GHz), φ, T₂ (ns), C]

% fit optimisation
opts = optimoptions('lsqcurvefit','Display','off');
% lsq curve fit is modifier to also return jacobian
[b_fit, ~, residual, ~,~,~, J]= lsqcurvefit(damped_sine, b0, tau, A_data, [], [], opts);

% estimate covariance matrix (for uncertainties)
mse = sum(residual.^2)/ (length(A_data)-length(b_fit));
covar = mse * inv(J' * J);

% standard errors i.e. 1 sigma uncertainties
% and also converting param std to full array before using
param_std = full(sqrt(diag(covar)));


% smooth fit curve
tau_fit = linspace(min(tau), max(tau), 1000);
A_fit = damped_sine(b_fit, tau_fit);

% plotting the fit
plot(tau_fit, A_fit, 'r-', 'LineWidth', 2);
legend('Data', 'Damped Sine Fit');

% Output fit params
fprintf('Fitted A0: %.3f ± %.3f \n', b_fit(1), param_std(1));
fprintf('Fitted Ω_R: %.6f ± %.6f GHz\n', b_fit(2), param_std(2));
fprintf('Fitted phase φ: %.3f ± %.3f rad\n', b_fit(3), param_std(3));
fprintf('Fitted T₂^{Rabi}: %.1f ± %.1f ns\n', b_fit(4), param_std(4));
fprintf('Fitted offset C: %.3f ± %.3f \n', b_fit(5), param_std(5));

