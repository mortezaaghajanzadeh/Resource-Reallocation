params = struct( ...
    'beta', 0.6, ...
    'gamma', 2.96, ...
    'alpha', 0.39, ...
    'A_hat', 0.1, ...
    'A_tilde', 0.03, ...
    'L', 250, ...
    'mu', 0, ...
    'theta', 1, ...
    'eta', Inf, ... % MATLAB uses 'Inf' for infinity
    'sigma', Inf, ... % Similarly, 'Inf' is used for infinity
    'r_b', 0.05, ...
    'green_premium', 0.0, ...
    'w', 500, ... % TSEK
    'tau_E', 0, ... % per ton
    'n', 1200 ... % Integer value for 'n'
);
% Calculate z_k using the defined parameters
z_k = function_z_k(params.alpha, params.gamma, params.r_b, params.green_premium, params.tau_E, params.A_tilde);

% Calculate z_l using the calculated z_k and other parameters
z_l = function_z_l(z_k, params.beta, params.alpha, params.gamma, params.w, params.green_premium, params.r_b);

% Display z_k and z_l
fprintf('z_k: %f\n', z_k);
fprintf('z_l: %f\n', z_l);

% Calculate and display intensity
intensity = function_intensity(params.A_tilde, params.A_hat, params.alpha, params.gamma, z_l, z_k, params.beta);
fprintf('intensity: (KG/SEK) %f\n', intensity);

% Calculate and display production
production = production_function(params.A_hat, params.beta, z_l, params.L) / 1000; % Converting to MSEK
fprintf('production: (MSEK) %f\n', production);

% Calculate and display emission
emission = intensity * production_function(params.A_hat, params.beta, z_l, params.L);
fprintf('emission: (ton) %f\n', emission);
fprintf('--------------------------------------------------\n');

% Simulation of multiple firms
[emission, production, intensity, ~, emission_cost] = simulate_firms(params.n, params.A_tilde, params.A_hat, params.L, params.alpha, params.gamma, params.r_b, params.green_premium, params.tau_E, params.beta, params.w, params.sigma);

% Display results from simulation
fprintf('Intensity: (KG/SEK) %f\n', intensity);
fprintf('Emission: (MTon) %f\n', emission / 1e6);
fprintf('Production: (BSek) %f\n', production / 1e6);
fprintf('Emission/Production: (KG/SEK) %f\n', emission / production);
fprintf('Emission cost: %f\n', sum(emission_cost));



function F = mySystem(x)
    simulate_firms(params.n, params.A_tilde, params.A_hat, params.L, params.alpha, params.gamma, params.r_b, params.green_premium, params.tau_E, params.beta, params.w, params.sigma)
end

function z_k = function_z_k(alpha, gamma, r_b, green_premium, tau_E, A_tilde)
    r_g = (1 - green_premium) * r_b;
    z_k = (alpha / (1 - alpha) * (r_b + tau_E / 1000 * A_tilde) / r_g) ^ (1 / gamma);
end
function z_l = function_z_l(z_k, beta, alpha, gamma, w, green_premium, r_b)
    r_g = (1 - green_premium) * r_b;
    z_l = (((1 - beta) / beta) / alpha) * ((alpha + (1 - alpha) * (z_k ^ (1 - gamma))) ^ (1 / (1 - gamma))) * (r_g / w);
end
function intensity = function_intensity(A_tilde, A_hat, alpha, gamma, z_l, z_k, beta)
    intensity = (A_tilde / A_hat) * (alpha * (z_k ^ (gamma - 1)) + (1 - alpha)) ^ (gamma / (1 - gamma)) * (z_l ^ (1 - beta));
end

function Y = production_function(A_hat, beta, z_l, L)
    Y = A_hat * (z_l ^ (-beta)) * L;
end
function [total_emissions, total_production, average_intensity, production, emission_cost] = simulate_firms(n, A_tilde, A_hat, L, alpha, gamma, r_b, green_premium, tau_E, beta, w, sigma)
    rng(0);  % Seed for reproducibility
    A_tilde_vector = (1 + lognrnd(0, 1, [n, 1])) * A_tilde;
    A_hat_vector = (1 + lognrnd(0, 1, [n, 1])) * A_hat;
    L_vector = lognrnd(-1.3, 1.8, [n, 1]) * L;
    intensity = zeros(n, 1);
    production = zeros(n, 1);
    emissions = zeros(n, 1);
    emission_cost = zeros(n, 1);

    for i = 1:n
        z_k = function_z_k(alpha, gamma, r_b, green_premium, tau_E, A_tilde_vector(i));
        z_l = function_z_l(z_k, beta, alpha, gamma, w, green_premium, r_b);
        IN = function_intensity(A_tilde_vector(i), A_hat_vector(i), alpha, gamma, z_l, z_k, beta);
        Y = production_function(A_hat_vector(i), beta, z_l, L_vector(i));
        intensity(i) = IN;
        production(i) = Y;
        emissions(i) = IN * Y;
        emission_cost(i) = tau_E / 1000 * IN * Y;
    end
    
    total_emissions = CES_aggregator(emissions, Inf);
    total_production = CES_aggregator(production, sigma);
    average_intensity = mean(intensity);
end

function result = CES_aggregator(array, sigma)
    if sigma ~= Inf
        result = sum(array .^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
    else
        result = sum(array);
    end
end