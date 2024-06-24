import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import statsmodels.api as sm



# set the F.O.C. functions

# $$ z_k = \frac{G}{B} = \left(\frac{\alpha}{1-\alpha}\frac{r^B}{r^G}\right)^{\gamma_s}$$
# $$ z_l = \frac{L}{K} =  \frac{1-\beta}{\beta} \frac{1}{\alpha}(\alpha_s + (1-\alpha_s) z_k^{-\frac{\gamma_s-1}{\gamma_s}})^{\frac{1}{\gamma_s-1}} \frac{r^G}{W} $$

def function_z_k(α,γ,r_b, green_premium,τ_E,A_tilde):
    r_g = (1-green_premium) * r_b
    return (α/(1-α) * (r_b +τ_E/1000*A_tilde )/r_g) ** (γ)

def function_green_ratio(α,γ,z_k):
    return (α + (1-α)*(z_k** (-(γ-1)/γ)) ) ** (-(γ/(γ-1)))

def function_brown_ratio(α,γ,z_k):
    return (α * (z_k ** ((γ-1)/γ)) + (1 - α)) ** (-(γ/(γ-1)))

def function_r(α,γ,z_k,green_premium,r_b):
    r_g = (1-green_premium) * r_b
    return function_green_ratio(α,γ,z_k) * r_g + function_brown_ratio(α,γ,z_k) * r_b

def function_z_l(z_k,β,α,γ,w,green_premium,r_b):
    r_g = (1-green_premium) * r_b
    return ((1-β)/β/α)* (function_green_ratio(α,γ,z_k)** (1/γ)) * (r_g/w)

def function_price_detail(A_tilde,A_hat,α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E):
    r_g = (1-green_premium) * r_b
    C_G = r_g * ( function_green_ratio(α,γ,z_k) * (z_l ** (β-1)))
    C_B = r_b * (function_brown_ratio(α,γ,z_k)) * (z_l ** (β-1))
    c_L = w * z_l ** (β)
    C_E = τ_E/1000 * A_tilde * C_B/r_b
    return C_G,C_B,c_L,C_E

def function_price(A_tilde,A_hat,α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E,detail=False):
    C_G,C_B,c_L,C_E = function_price_detail(A_tilde,A_hat,α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E)
    if detail:
        return C_G,C_B,c_L,C_E, σ/(σ-1) * (C_G + C_B + c_L + C_E) / A_hat
    else:
        return σ/(σ-1) * (C_G + C_B + c_L + C_E) / A_hat

def function_intensity(A_tilde,A_hat,α,γ,z_l,z_k,β): # This the E/Y and different from E/PY
    return (A_tilde/A_hat) * function_brown_ratio(α,γ,z_k) * (z_l ** (β-1))

    
def optimal_labor(A_tilde,A_hat,α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E,κ=0.1):
    p = function_price(A_tilde,A_hat,α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E)
    return κ * z_l ** (β) * (p ** (-σ))/A_hat

def production_function(A_hat,β,z_l,L):
    return A_hat * (z_l ** (β)) * L

def brown_capital(A_hat,α,z_k,z_l,γ,β,Y):
    return (Y/A_hat) * ((α * (z_k ** ((γ-1)/γ)) + (1 - α)) ** (γ / (1-γ))) * (z_l ** (β-1))

def green_capital(A_hat,α,z_k,z_l,γ,β,Y):
    return (Y/A_hat) * ((α + (1-α)*(z_k** (-(γ-1)/γ)) ) ** (γ / (1-γ))) * (z_l ** (β-1))


def optimal_ratios(input_0,τ_E):
    z_k = function_z_k(input_0['α'],input_0['γ'],input_0['r_b'], input_0['green_premium'],τ_E,input_0['A_tilde'])
    z_l = function_z_l(z_k,input_0['β'],input_0['α'],input_0['γ'],input_0['w'],input_0['green_premium'],input_0['r_b'])
    return z_k,z_l

def ratios_gen(input_0):
    z_k,z_l = optimal_ratios(input_0,input_0['τ_E'])
    l = optimal_labor(input_0['A_tilde'],input_0['A_hat'],input_0['α'],input_0['γ'],z_l,z_k,input_0['β'],input_0['w'],input_0['green_premium'],input_0['r_b'],input_0['σ'],input_0['τ_E'])
    y = production_function(input_0['A_hat'],input_0['β'],z_l,l)
    p = function_price(input_0['A_tilde'],input_0['A_hat'],input_0['α'],input_0['γ'],z_l,z_k,input_0['β'],input_0['w'],input_0['green_premium'],input_0['r_b'],input_0['σ'],input_0['τ_E'])
    b = brown_capital(input_0['A_hat'],input_0['α'],z_k,z_l,input_0['γ'],input_0['β'],y)
    g = green_capital(input_0['A_hat'],input_0['α'],z_k,z_l,input_0['γ'],input_0['β'],y)
    return z_k,z_l,l,y*1e3,p, g, b




def simulate_firms(par):
    α = par['α']
    β = par['β']
    γ = par['γ']
    σ = par['σ']
    τ_E = par['τ_E']
    green_premium = par['green_premium']
    r_b = par['r_b']
    w = par['w']
    A_tilde = par['A_tilde']
    A_hat = par['A_hat']
    n = par['n']
    sd_hat = par['sd_hat']
    sd_tilde = par['sd_tilde']
    rho = par['rho']
    cov = rho * sd_hat * sd_tilde

    # simulate over multiple firms
    r_g = (1-green_premium) * r_b
    np.random.seed(0)
    # A_tilde_vector = (1 + np.random.lognormal(mean=0,sigma=sd,size=n)) * A_tilde
    # A_hat_vector = (1 + np.random.lognormal(mean=0,sigma=sd,size=n)) * A_hat

    A_tilde_mu = np.log(A_tilde)
    A_hat_mu = np.log(A_hat)
    #build two correlated vectors
    mu = np.array([A_tilde_mu,A_hat_mu])

    # The desired covariance matrix.
    r = np.array([
            [sd_tilde ** 2, cov],
            [cov, sd_hat ** 2]
        ])
    # np.corrcoef(y[:,0],y[:,1])
    # Generate the random samples.
    rng = np.random.default_rng(seed=0)
    A_vector = rng.multivariate_normal(mu, r, size=n)
    A_vector = np.exp(A_vector)


    intensity = []
    production = []
    emissions = []
    emission_cost = []
    G_c = []
    B_c = []
    labor = []
    income = []
    cost_share = []
    price = []
    z_l = []
    z_k = []
    K = []
    parameters = par.copy()
    for i in A_vector:
        parameters['A_tilde'] = i[0]
        parameters['A_hat'] = i[1]
        k,l_ratio,l,Y,p,g, b = ratios_gen(parameters)
        z_l.append(l_ratio)
        z_k.append(k)
        emission = i[0] * b
        IN = emission / p / Y
        intensity.append(IN)
        production.append(Y)
        emissions.append(emission)
        emission_cost.append(τ_E/1000 * emission)
        G_c.append(g)
        B_c.append(b)
        labor.append(l)
        income.append(p*Y)
        cost_share.append(τ_E * emission / p/Y)
        price.append(p)
        K.append(l/l_ratio)
    # result = (
    #     CES_aggregator(emissions,np.inf),
    #     CES_aggregator(production,σ),
    #     sum(intensity)/n,
    #     production,
    #     emission_cost,
    #     G_c,
    #     B_c,
    #     labor,
    #     income,
    #     cost_share,
    #     price
    #     )
    result = (
        np.array(emissions),
        np.array(production),
        np.array(intensity),
        np.array(emission_cost),
        np.array(G_c),
        np.array(B_c),
        np.array(labor),
        np.array(income),
        np.array(cost_share),
        np.array(price),
        np.array(z_k),
        np.array(z_l),
        np.array(K)
    )
    return result

def CES_aggregator(array,σ):
    if σ != np.inf:
        return sum([i ** ((σ-1)/σ) for i in array]) ** (σ/(σ-1))
    else:
        return sum(array)

