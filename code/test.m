clear
parameters = struct(...
    'beta', 0.6, ...
    'gamma', 2.7268, ...
    'alpha', 0.25, ...
    'A_hat', 2, ...
    'A_tilde', .0002, ...
    'mu', 0, ...
    'theta', 1, ...
    'eta', inf, ...
    'sigma',5, ...
    'r_b', 0.05, ...
    'green_premium', 0, ...
    'w', 0.500, ... % MSEK
    'tau_E', 0, ... % per ton
    'n', 1000, ...
    'sd_hat', 0.45, ...
    'sd_tilde', 2, ...
    'rho', 0.9 ...
);
clc
obj = functionsContainer();

% Define tax values
tax_1 = 0;
tax_2 = 100;
tax_3 = 1300;

% Generate tables for each tax scenario
df1 = obj.gen_df(parameters, tax_1);
df1.Properties.RowNames = strcat('τ_E = ', string(tax_1), '_', df1.Properties.RowNames); % Add prefix to row names

df2 = obj.gen_df(parameters, tax_2);
df2.Properties.RowNames = strcat('τ_E = ', string(tax_2), '_', df2.Properties.RowNames); % Add prefix to row names

df3 = obj.gen_df(parameters, tax_3);
df3.Properties.RowNames = strcat('τ_E = ', string(tax_3), '_', df3.Properties.RowNames); % Add prefix to row names

% Concatenate tables
df = [df1; df2; df3]
