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

def function_intensity(A_tilde,A_hat,α,γ,z_l,z_k,β):
    return (A_tilde/A_hat) * function_brown_ratio(α,γ,z_k) * (z_l ** (β-1))

    
def optimal_labor(A_tilde,A_hat,α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E,κ=1e3):
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




def simulate_firms(n,A_tilde,A_hat,α,γ,r_b, green_premium,τ_E,β,w,σ):
    # simulate over multiple firms
    r_g = (1-green_premium) * r_b
    np.random.seed(0)
    A_tilde_vector = (1 + np.random.lognormal(mean=0,sigma=1,size=n)) * A_tilde
    A_hat_vector = (1 + np.random.lognormal(mean=0,sigma=1,size=n)) * A_hat
    intensity = []
    production = []
    emissions = []
    emission_cost = []
    G_c = []
    B_c = []
    labor = []
    for i in zip(A_tilde_vector,A_hat_vector):
        z_k = function_z_k(α,γ,r_b, green_premium,τ_E,i[0])
        z_l = function_z_l(z_k,β,α,γ,w,green_premium,r_b)
        p = function_price(i[0],i[1],α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E)
        l = optimal_labor(i[0],i[1],α,γ,z_l,z_k,β,w,green_premium,r_b,σ,τ_E)
        IN = function_intensity(i[0],i[1],α,γ,z_l,z_k,β)
        Y = production_function(i[1],β,z_l,l)
        intensity.append(IN)
        production.append(Y)
        emissions.append(IN * Y)
        emission_cost.append(τ_E/1000 * IN * Y)
        G_c.append(green_capital(i[1],α,z_k,z_l,γ,β,Y))
        B_c.append(brown_capital(A_hat,α,z_k,z_l,γ,β,Y))
        labor.append(l)
    return CES_aggregator(emissions,np.inf),CES_aggregator(production,σ),sum(intensity)/n,production,emission_cost,G_c,B_c,labor
def CES_aggregator(array,σ):
    if σ != np.inf:
        return sum([i ** ((σ-1)/σ) for i in array]) ** (σ/(σ-1))
    else:
        return sum(array)

def gen_df(params,τ_E):
    params_0 = params.copy()
    params_0['τ_E'] = τ_E 
    params_1 = params_0.copy()
    res_1 = ratios_gen(params_1)
    params_2 = params_0.copy()
    params_2['A_hat'] = params_2['A_hat'] * 0.8
    params_2['A_tilde'] = params_2['A_tilde'] * 1.2
    res_2 = ratios_gen(params_2)
    # put the results into a dataframe
    df = pd.DataFrame([res_1,res_2],columns=['z_k','z_l','l','y','p','g','b'])
    df['income'] = (df.y * df.p).round(2)
    # df['emission'] = df.emission.round(2)
    # add sum row
    df.loc['sum'] = np.nan
    df.loc['sum','l']= df.l.sum()
    df.loc['sum','g']= df.g.sum()
    df.loc['sum','b']= df.b.sum()
    df.loc['sum','income']= df.income.sum()
    p_s = 0
    y_s = 0
    for i in [res_1,res_2]:
        p_s += i[4] ** (1-params_0['σ'])
        y_s += i[3] ** ((params_0['σ']-1)/(params_0['σ']))
    p_s = p_s ** (1/(1-params_0['σ']))
    y_s = y_s ** (params_0['σ']/(params_0['σ']-1))
    df.loc['sum','p'] = p_s
    # df.p = df.p/p_s
    df.loc['sum','y'] = y_s
    return df
