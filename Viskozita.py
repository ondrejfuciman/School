###########################################
# Parametry
###########################################

#průměr kuliček
d = 2/1000 #m
sigma_d = 0.005/1000 #

#hmotnost (ložiskových) kuliček
m = ((0.2959/9 + 0.7824/24)/2)/1000 #kg

#odchylka hmotnosti
sigma_m = (4*10**(-4))/1000 #kg

#hustota ložiskových kuliček podle dokumentu z praktik
rho_k = 7860 #kg/m**3

#tíhové zrychlení - přesnější než čistě 9,81 a mezinárodně uznávaná
import scipy
g = scipy.constants.g #m*s**(-2)

#vertikální vzdálenost padání kuliček v glycerinu
h_g = 53.0/100 #m 

#časy dopadu kuliček v glycerinu
t_g = [
    31.39,
    31.45,
    32.56,
    31.10,
    31.41,
    32.28,
    31.47,
    31.79,
    31.16,
    31.68,
    31.68,
    31.31,
    31.71,
    31.82
]

#hustota glycerinu
rho_g = 1230 # kg/m**3

#chyba odhadu hustoty glycerinu
sigma_rho_g = 10 #kg/m**3

#vertikální vzdálenost padání kuliček v oleji
h_o = 49.7/100 #m 

#časy dopadu kuliček v oleji
t_o = [
    29.35,
    29.40,
    29.32,
    29.36,
    29.28,
    29.30,
    29.36,
    29.43,
    29.35,
    29.38,
    29.51,
    29.46
]

sigma_t = 0.5/1000 #s

#hustota oleje
rho_o = 955 #kg/m**3

#chyba odhadu hustoty oleje
sigma_rho_o = 5 #kg/m**3

#reakční faktor - Ondra
sigma_t_r_O = 0.037 #s

#reakční faktor - Artem
sigma_t_r_A = 0.018 #s

#nejistota vertikální vzdálenosti padání
sigma_h = 0.05/100 #m

import numpy as np
import sympy as sp

def dynamická_viskozita(rho_l:float,  
    sigma_rho_l:float, t:list, h:float, sigma_t:float)->tuple[float, float]:

    t_avg = np.mean(t)
    sigma_t_stat = np.std(t, ddof=1)

    sigma_t = np.sqrt(sigma_t_stat**2 + sigma_t**2)

    d_s, g_s, rho_k_s, rho_l_s, t_avg_s, h_s = sp.symbols('d g rho_k rho_l t h')

    eta_expr = (2/9)*g_s*(d_s/2)**2*(rho_k_s-rho_l_s)*(t_avg_s/h_s)

    collapsed_to_value = {
        d_s: d,
        g_s: g,
        rho_k_s: rho_k,
        rho_l_s: rho_l,
        t_avg_s: t_avg,
        h_s: h,
    }

    eta_value = float(eta_expr.subs(collapsed_to_value))


    partial_eta_partial_d = float(sp.diff(eta_expr, d_s).subs(collapsed_to_value))
    partial_eta_partial_h = float(sp.diff(eta_expr, h_s).subs(collapsed_to_value))
    partial_eta_partial_rho_l = float(sp.diff(eta_expr, rho_l_s).subs(collapsed_to_value))
    partial_eta_partial_t_avg = float(sp.diff(eta_expr, t_avg_s).subs(collapsed_to_value))

    sigma_eta = np.sqrt((partial_eta_partial_d*sigma_d)**2 +
                        (partial_eta_partial_h*sigma_h)**2 +
                        (partial_eta_partial_rho_l*sigma_rho_l)**2 +
                        (partial_eta_partial_t_avg*sigma_t)**2)
    
    return eta_value, sigma_eta

eta_g1, sigma_eta_g1 = dynamická_viskozita(rho_g, sigma_rho_g, t_g, h_g, sigma_t_r_O)
eta_o1, sigma_eta_o1 = dynamická_viskozita(rho_o, sigma_rho_o, t_o, h_o, sigma_t_r_A)

print(r"Dynamická viskozita glycerinu s $\rho_k$ opsaným z dokumnetu vychází", eta_g1,"±", sigma_eta_g1)
print(r"Dynamická viskozita ricinového oleje s $\rho_k$ opsaným z dokumnetu vychází", eta_o1,"±", sigma_eta_o1)

#Vážení kuliček by bylo zbytečné, kdybychom si vzali hustotu čistě z dokumentu

def dynamická_viskozita_měřená_rho_k(rho_l:float,  
    sigma_rho_l:float, t:list, h:float, sigma_t:float)->tuple[float, float]:

    t_avg = np.mean(t)
    sigma_t_stat = np.std(t, ddof=1)

    C = (4/3 * np.pi * (1/2)**3)**(-1)

    sigma_t = np.sqrt(sigma_t_stat**2 + sigma_t**2)

    C_s, d_s, g_s, m_s, rho_l_s, t_avg_s, h_s = sp.symbols('C d g m rho_l t h')

    eta_expr = (2/9)*g_s*(d_s/2)**2*(C*m_s*d_s**(-3)-rho_l_s)*(t_avg_s/h_s)

    collapsed_to_value = {
        C_s: C,
        d_s: d,
        g_s: g,
        m_s: m,
        rho_l_s: rho_l,
        t_avg_s: t_avg,
        h_s: h,
    }

    eta_value = float(eta_expr.subs(collapsed_to_value))

    partial_eta_partial_d = float(sp.diff(eta_expr, d_s).subs(collapsed_to_value))
    partial_eta_partial_h = float(sp.diff(eta_expr, h_s).subs(collapsed_to_value))
    partial_eta_partial_m = float(sp.diff(eta_expr, m_s).subs(collapsed_to_value))
    partial_eta_partial_rho_l = float(sp.diff(eta_expr, rho_l_s).subs(collapsed_to_value))
    partial_eta_partial_t_avg = float(sp.diff(eta_expr, t_avg_s).subs(collapsed_to_value))
    

    sigma_eta = np.sqrt((partial_eta_partial_d*sigma_d)**2 +
                        (partial_eta_partial_h*sigma_h)**2 +
                        (partial_eta_partial_m*sigma_m)**2 +
                        (partial_eta_partial_rho_l*sigma_rho_l)**2 +
                        (partial_eta_partial_t_avg*sigma_t)**2)
    
    return eta_value, sigma_eta

eta_g2, sigma_eta_g2 = dynamická_viskozita_měřená_rho_k(rho_g, sigma_rho_g, t_g, h_g, sigma_t_r_O)
eta_o2, sigma_eta_o2 = dynamická_viskozita_měřená_rho_k(rho_o, sigma_rho_o, t_o, h_o, sigma_t_r_A)

print("Dynamická viskozita glycerinu s dopočteným $\rho_k$ vychází", eta_g2,"±", sigma_eta_g2)
print("Dynamická viskozita ricinového oleje s dopočteným $\rho_k$ vychází", eta_o2,"±", sigma_eta_o2)

print(m)

def vlastní_hustota():
    C_s, m_s, d_s = sp.symbols('C m d')

    C = (4/3 * np.pi * (1/2)**3)**(-1)

    rho_expr = C_s * m_s * d_s**(-3)

    vals = {
        C_s: C,
        m_s: m,
        d_s: d
    }

    rho_val = float(rho_expr.subs(vals))

    partial_m = float(sp.diff(rho_expr, m_s).subs(vals))
    partial_d = float(sp.diff(rho_expr, d_s).subs(vals))

    sigma_rho = np.sqrt(
        (partial_m * sigma_m)**2 +
        (partial_d * sigma_d)**2
    )

    return rho_val, sigma_rho

rho_k_calc, sigma_rho_calc = vlastní_hustota()

print("Vlastní dopočtená hustota vyšla", rho_k_calc,"±", sigma_rho_calc)


print(g)

sigma_t_g_stat = np.std(t_g, ddof=1)
sigma_t_g = np.sqrt(sigma_t_g_stat**2 + sigma_t_r_O**2)

sigma_t_o_stat = np.std(t_o, ddof=1)
sigma_t_o = np.sqrt(sigma_t_o_stat**2 + sigma_t_r_A**2)

print("Chyba doby padání v glycerolu je", sigma_t_g, "a v oleji je ", sigma_t_o)
