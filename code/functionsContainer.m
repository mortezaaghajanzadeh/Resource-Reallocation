classdef functionsContainer
    methods
        function [z_k] = function_z_k(obj, alpha, gamma, r_b, green_premium, tau_E, A_tilde)
            r_g = (1 - green_premium) * r_b;
            z_k = (alpha / (1 - alpha) * (r_b + tau_E * A_tilde) / r_g) .^ gamma;
        end

        function [ratio] = function_green_ratio(obj, alpha, gamma, z_k) % G/K
            ratio = (alpha + (1 - alpha) * (z_k .^ (-(gamma - 1) / gamma))) .^ (-(gamma / (gamma - 1)));
        end

        function [ratio] = function_brown_ratio(obj, alpha, gamma, z_k) % B/K
            ratio = (alpha * (z_k .^ ((gamma - 1) / gamma)) + (1 - alpha)) .^ (-(gamma / (gamma - 1)));
        end

        function [r] = function_r(obj, alpha, gamma, z_k, green_premium, r_b)
            r_g = (1 - green_premium) * r_b;
            r = obj.function_green_ratio(alpha, gamma, z_k) * r_g + obj.function_brown_ratio(alpha, gamma, z_k) * r_b;
        end

        function [z_l] = function_z_l(obj, z_k, beta, alpha, gamma, w, green_premium, r_b) % L/K
            r_g = (1 - green_premium) * r_b;
            z_l = ((1 - beta) / beta / alpha) * (obj.function_green_ratio(alpha, gamma, z_k) .^ (1 / gamma)) * (r_g / w);
        end

        function [C_G, C_B, c_L, C_E] = function_price_detail(obj, A_tilde, A_hat, alpha, gamma, z_l, z_k, beta, w, green_premium, r_b, sigma, tau_E)
            r_g = (1 - green_premium) * r_b;
            C_G = r_g * (obj.function_green_ratio(alpha, gamma, z_k) * (z_l .^ (beta - 1)));
            C_B = r_b * (obj.function_brown_ratio(alpha, gamma, z_k)) * (z_l .^ (beta - 1));
            c_L = w * z_l .^ beta;
            C_E = tau_E * A_tilde * C_B / r_b;
        end

        function [result] = function_price(obj, A_tilde, A_hat, alpha, gamma, z_l, z_k, beta, w, green_premium, r_b, sigma, tau_E, detail)
            [C_G, C_B, c_L, C_E] = obj.function_price_detail(A_tilde, A_hat, alpha, gamma, z_l, z_k, beta, w, green_premium, r_b, sigma, tau_E);
            if detail
                result = [C_G, C_B, c_L, C_E, sigma / (sigma - 1) * (C_G + C_B + c_L + C_E) / A_hat];
            else
                result = sigma / (sigma - 1) * (C_G + C_B + c_L + C_E) / A_hat;
            end
        end

        function [intensity] = function_intensity(obj, A_tilde, A_hat, alpha, gamma, z_l, z_k, beta) % This is E/Y and different from E/PY
            intensity = (A_tilde / A_hat) * obj.function_brown_ratio(alpha, gamma, z_k) * (z_l .^ (beta - 1));
        end

        function [labor] = optimal_labor(obj, A_tilde, A_hat, alpha, gamma, z_l, z_k, beta, w, green_premium, r_b, sigma, tau_E, kappa)
            p = obj.function_price(A_tilde, A_hat, alpha, gamma, z_l, z_k, beta, w, green_premium, r_b, sigma, tau_E, false);
            labor = kappa * z_l .^ beta * (p .^ (-sigma)) / A_hat;
        end

        function [Y] = production_function(obj, A_hat, beta, z_l, L)
            Y = A_hat * (z_l .^ (-beta)) * L;
        end

        function [K_B] = brown_capital(obj, A_hat, alpha, z_k, z_l, gamma, beta, Y)
            % K_B = (Y / A_hat) * ((alpha * (z_k .^ ((gamma - 1) / gamma)) + (1 - alpha)) .^ (gamma / (1 - gamma))) * (z_l .^ (beta - 1));
            K_B = (Y / A_hat) * obj.function_brown_ratio(alpha, gamma, z_k) * (z_l .^ (beta - 1));
        end

        function [K_G] = green_capital(obj, A_hat, alpha, z_k, z_l, gamma, beta, Y)
            % K_G = (Y / A_hat) * ((alpha + (1 - alpha) * (z_k .^ (-(gamma - 1) / gamma))) .^ (gamma / (1 - gamma))) * (z_l .^ (beta - 1));
            K_G = (Y / A_hat) * obj.function_green_ratio(alpha, gamma, z_k) * (z_l .^ (beta - 1));
        end

        function [z_k, z_l] = optimal_ratios(obj, input_0, tau_E)
            z_k = obj.function_z_k(input_0.alpha, input_0.gamma, input_0.r_b, input_0.green_premium, tau_E, input_0.A_tilde);
            z_l = obj.function_z_l(z_k, input_0.beta, input_0.alpha, input_0.gamma, input_0.w, input_0.green_premium, input_0.r_b);
        end

        function [z_k, z_l, l, Y, p, g, b] = ratios_gen(obj, input_0)
            [z_k, z_l] = obj.optimal_ratios(input_0, input_0.tau_E);
            l = obj.optimal_labor(input_0.A_tilde, input_0.A_hat, input_0.alpha, input_0.gamma, z_l, z_k, input_0.beta, input_0.w, input_0.green_premium, input_0.r_b, input_0.sigma, input_0.tau_E, 0.1);
            Y = obj.production_function(input_0.A_hat, input_0.beta, z_l, l);
            p = obj.function_price(input_0.A_tilde, input_0.A_hat, input_0.alpha, input_0.gamma, z_l, z_k, input_0.beta, input_0.w, input_0.green_premium, input_0.r_b, input_0.sigma, input_0.tau_E, false);
            b = obj.brown_capital(input_0.A_hat, input_0.alpha, z_k, z_l, input_0.gamma, input_0.beta, Y);
            g = obj.green_capital(input_0.A_hat, input_0.alpha, z_k, z_l, input_0.gamma, input_0.beta, Y);
        end
        
        function result = CES_aggregator(obj,array, sigma)
            if sigma ~= inf
                result = sum(array .^ ((sigma - 1) / sigma)) .^ (sigma / (sigma - 1));
            else
                result = sum(array);
            end
        end

        function result = simulate_firms(obj,par)
            % alpha = par.alpha;
            % beta = par.beta;
            % gamma = par.gamma;
            % sigma = par.sigma;
            tau_E = par.tau_E;
            % green_premium = par.green_premium;
            % r_b = par.r_b;
            % w = par.w;
            A_tilde = par.A_tilde;
            A_hat = par.A_hat;
            n = par.n;
            sd_hat = par.sd_hat;
            sd_tilde = par.sd_tilde;
            rho = par.rho;
            cov = rho * sd_hat * sd_tilde;
        
            % Generate correlated log-normal distributions
            mu = [log(A_tilde), log(A_hat)];
            R = [sd_tilde^2, cov; cov, sd_hat^2];
            rng(0); % Seed for reproducibility
            A_vector = mvnrnd(mu, R, n);
            A_vector = exp(A_vector); % Convert to lognormal
        
            % Initialize arrays
            intensity = zeros(1, n);
            production = zeros(1, n);
            emissions = zeros(1, n);
            emission_cost = zeros(1, n);
            G_c = zeros(1, n);
            B_c = zeros(1, n);
            labor = zeros(1, n);
            income = zeros(1, n);
            cost_share = zeros(1, n);
            price = zeros(1, n);
            z_l = zeros(1, n);
            z_k = zeros(1, n);
            K = zeros(1, n);
        
            % Simulate for each firm
            for i = 1:n
                parameters = par; % Copy parameters
                parameters.A_tilde = A_vector(i, 1);
                parameters.A_hat = A_vector(i, 2);
                [k, l_ratio, l, Y, p, g, b] = ratios_gen(obj, parameters);
                z_l(i) = l_ratio;
                z_k(i) = k;
                emission = A_vector(i, 1) * b;
                IN = emission / (p * Y);
                intensity(i) = IN;
                production(i) = Y;
                emissions(i) = emission;
                emission_cost(i) = tau_E * emission;
                G_c(i) = g;
                B_c(i) = b;
                labor(i) = l;
                income(i) = p * Y;
                cost_share(i) = tau_E * emission / (p * Y);
                price(i) = p;
                K(i) = l / l_ratio;
            end
        
            result = struct('emissions', emissions, 'production', production, 'intensity', intensity, ...
                            'emission_cost', emission_cost, 'G_c', G_c, 'B_c', B_c, 'labor', labor, ...
                            'income', income, 'cost_share', cost_share, 'price', price, 'z_k', z_k, ...
                            'z_l', z_l, 'K', K);
        end

        function df = gen_df(obj,params, tau_E)
            columns =  {'z_k', 'z_l', 'l', 'y', 'p', 'k', 'g', 'b', 'e', 'coef'};
            params_1 = params;
            params_1.tau_E = tau_E;
            [z_k, z_l, l, y, p, g, b] = ratios_gen(obj,params_1);
            e = b * params_1.A_tilde; % Assuming 'l * params_1.A_tilde' is intended to compute 'e'
            coef = function_coef(obj,params_1, tau_E);
            k = l / z_l;
            
            % Store the first row of results
            row1 = table(z_k, z_l, l, y, p, k, g, b, e, coef, 'VariableNames',columns);
            
            % Copy the parameters and modify for the second scenario
            params_2 = params;
            params_2.tau_E = tau_E;
            params_2.A_hat = params_2.A_hat * 0.8;
            params_2.A_tilde = params_2.A_tilde * 1.2;
            [z_k, z_l, l, y, p, g, b] = ratios_gen(obj,params_2);
            e = b * params_2.A_tilde;
            coef = function_coef(obj,params_2, tau_E);
            k = l / z_l;
            
            % Store the second row of results
            row2 = table(z_k, z_l, l, y, p, k, g, b, e, coef, 'VariableNames', columns);
            
            % Create a MATLAB table from the results
            df = [row1; row2];
            % Calculate derived columns and round 'income'
            df.income = round(df.y .* df.p, 2);
            df.intensity = (df.e ./ (df.y .* df.p));
            
            % Add a sum row
            l_s = sum(df.l);
            k_s = sum(df.k);
            g_s = sum(df.g);
            b_s = sum(df.b);
            e_s = sum(df.e);
            income_s = sum(df.income);
            p_s = sum(df.p .^ (1 - params.sigma));
            y_s = sum(df.y .^ ((params.sigma - 1) / params.sigma));
            p_s = p_s ^ (1 / (1 - params.sigma));
            y_s = y_s ^ (params.sigma / (params.sigma - 1));

            %Append sum results to the table
            df = [df; array2table([nan nan l_s y_s p_s k_s g_s b_s e_s nan income_s nan], ...
                'VariableNames', [columns {'income','intensity'}])];
            df.Properties.RowNames(1) = {'1'};
            df.Properties.RowNames(2) = {'2'};
            df.Properties.RowNames(3) = {'sum'};
        end

        function coef = function_coef(obj,parameters, tau_E)
            A_tilde = parameters.A_tilde;
            alpha = parameters.alpha;
            gamma = parameters.gamma;
            [z_l, z_k] = optimal_ratios(obj,parameters, tau_E);
            green_premium = parameters.green_premium;
            r_b = parameters.r_b;
            w = parameters.w;
            beta = parameters.beta;
            r = function_r(obj,alpha, gamma, z_k, green_premium, r_b);

            if tau_E == 0
                coef = 0;
            else
                epsilon_change = A_tilde / (r_b + tau_E / 1000 * A_tilde) * (...
                    -gamma * alpha / (function_green_ratio(obj,alpha, gamma, z_k) ^ (-((gamma - 1) / gamma))) ...
                    - (1 - alpha) / (function_green_ratio(obj,alpha, gamma, z_k) ^ (-((gamma - 1) / gamma))) * ...
                    ( 2 * beta - (1 +  r/z_l/(w + r/z_l))));
                c_change = 1 / (tau_E / 1000) + epsilon_change;
                coef = epsilon_change / c_change;
            end
        end

        function K = function_K(obj,g,b,alpha,gamma)
            K = (alpha * g ^ ((gamma - 1) / gamma) + (1 - alpha) * b ^ ((gamma - 1) / gamma)) ^ (gamma / (gamma - 1));
        end

        function res_table = simulate_firms_table(obj,parameters)
            res = obj.simulate_firms(parameters);
            % Extract results from the structure
            emissions = res.emissions;
            production = res.production;
            intensity = res.intensity;
            emission_cost = res.emission_cost;
            G_c = res.G_c;
            B_c = res.B_c;
            l = res.labor;
            income = res.income;
            cost_share = res.cost_share;
            price = res.price;
            z_k_ratios = res.z_k;
            z_l_ratios = res.z_l;
            K = res.K;
            wage_share = (res.labor * parameters.w) ./ res.income;

            % Aggregating with CES_aggregator function
            Y = obj.CES_aggregator(production, parameters.sigma);
            e = sum(emissions);

            % Generating statistics and storing them in a MATLAB table
            res_table = table('Size', [14, 4], ...
                            'VariableTypes', repmat({'double'}, 1, 4), ...
                            'VariableNames', {'mean', 'median', 'std', 'sum'}, ...
                            'RowNames', {'emissions', 'production', 'intensity', 'emission_cost', 'G_c', 'B_c', 'K', 'l', 'income', 'cost_share', 'wage_share', 'price', 'z_k', 'z_l'});
            data_variables = {emissions, production, intensity, emission_cost, G_c, B_c, K, l, income, cost_share, wage_share, price, z_k_ratios, z_l_ratios};
            stat_functions = {@mean, @median, @std, @sum};
            stat_names = {'mean', 'median', 'std', 'sum'};

            for i = 1:numel(data_variables)
                for j = 1:numel(stat_functions)
                    if strcmp(stat_names{j}, 'sum') && any(strcmp(res_table.Properties.RowNames{i}, {'intensity', 'price', 'z_k', 'z_l'}))
                        % Handle cases where sum is not applicable
                        res_table.(stat_names{j})(i) = NaN;
                    else
                        res_table.(stat_names{j})(i) = stat_functions{j}(data_variables{i});
                    end
                end
            end

            % For specific cases where the sum needs to be replaced with special calculations
            res_table.sum('production') = Y;
            res_table.sum('emissions') = e;
        end
    end
end
