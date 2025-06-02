import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from scipy.stats import norm, multivariate_normal

# Constantes físicas
G = 6.67430e-8  # cm^3 g^-1 s^-2
Msun = 1.989e33  # g
Rsun = 6.957e10  # cm

# Parámetros estelares de TW Hya
M_star = 0.8 * Msun
R_star = 1.1 * Rsun
AU = 1.496e13 #cm
#R_trunc = 4 * AU  # radio de truncamiento típico
R_trunc = 3.5 * R_star  # radio de magnetosferico típico

def get_inv_cdf_mass_function(n=1000, m_min=3e-14, m_max=4.5e-11, alpha=-1.2):
    x = np.logspace(np.log10(m_min), np.log10(m_max), n)
    pdf = x**alpha
    pdf /= np.trapz(pdf, x)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return interp.interp1d(cdf, x, bounds_error=False, fill_value=(x[0], x[-1]))

# Nueva función: muestreo copulado en espacio log-log
def sample_loglog_copula_correlated(n_samples, logmass_range, logdur_range, rho=0.99):
    mean = [0, 0]
    cov = [[1, rho], [rho, 1]]
    mvn_samples = multivariate_normal.rvs(mean=mean, cov=cov, size=n_samples)
    u = norm.cdf(mvn_samples)
    # Volver a log10 escala física
    log_m_min, log_m_max = logmass_range
    log_t_min, log_t_max = logdur_range
    log_masses = log_m_min + u[:, 0] * (log_m_max - log_m_min)
    log_durations = log_t_min + u[:, 1] * (log_t_max - log_t_min)
    masses = 10**log_masses
    durations = 10**log_durations
    return masses, durations

def sample_masses_and_durations(method='gaussian', **kwargs):
    if method == 'gaussian':
        keys = ['n_samples', 'logmass_range', 'logdur_range', 'rho']
        filtered_kwargs = {k: kwargs[k] for k in keys if k in kwargs}
        return sample_loglog_copula_correlated(**filtered_kwargs)
    else:
        raise ValueError("Método no reconocido: usa 'gaussian'")

#def run_clumpy_model(alpha=-1.2, mu_log=np.log(1.5), sigma_log=1, rho=0.98, total_days=100, cadence_min=2):
def run_clumpy_model(alpha=-1.2, mu_log=np.log(1.5), sigma_log=1, rho=0.99, total_days=100, cadence_min=2):
    # Parámetros de simulación
    dt = cadence_min * 60  # segundos
    time = np.arange(0, total_days * 86400, dt)  # en segundos

    # Número esperado de bursts por campaña (~40-50)
    n_bursts = np.int64(47 * (total_days/100) * 4)

    inv_cdf_mass = get_inv_cdf_mass_function(alpha=alpha)
    # Crear distribución personalizada para duraciones
    x_dur = np.logspace(np.log10(0.045), np.log10(17), 1000)
    pdf_dur = (1 / x_dur) * np.exp(- (np.log(x_dur) - mu_log)**2 / (2 * sigma_log**2))
    pdf_dur /= np.trapz(pdf_dur, x_dur)
    cdf_dur = np.cumsum(pdf_dur)
    cdf_dur /= cdf_dur[-1]
    inv_cdf_dur = interp.interp1d(cdf_dur, x_dur, bounds_error=False, fill_value=(x_dur[0], x_dur[-1]))

    # Definir rangos en log10 para masas y duraciones
    logmass_range = (np.log10(3e-14), np.log10(4.5e-11))  # rango de masas en M_sun
    logdur_range = (np.log10(0.045), np.log10(17))        # rango de duración en días

    # Muestreo copulado en espacio log-log
    copula_method = 'gaussian'  # cambiar a 'gaussian' si se desea

    clump_masses, durations_days = sample_masses_and_durations(
        method=copula_method,
        n_samples=n_bursts,
        logmass_range=logmass_range,
        logdur_range=logdur_range,
        rho=rho  # solo aplica para Gaussiana
    )

    # Simulación de bursts
    lc = np.zeros_like(time, dtype=np.float64)
    burst_data = []

    rho_clump = 5e-12  # g/cm^3
    v_ff2 = np.sqrt(2 * G * M_star * (1/R_star - 1/R_trunc))
    durations_seconds2 = ((3 * clump_masses * Msun) / (4 * np.pi * rho_clump))**(1/3) / v_ff2
    durations_days2 = durations_seconds2 / 86400

    for i in range(n_bursts):
        t0 = np.random.uniform(time[0], time[-1])  # centro del burst
        m_clump = clump_masses[i] * Msun  # en gramos

        # velocidad de caída libre
        v_ff = np.sqrt(2 * G * M_star * (1/R_star - 1/R_trunc))

        # energía liberada
        E_acc = 0.5 * m_clump * v_ff**2  # en erg

        # duración típica: log-uniforme entre 0.5 y 3 días
        tau = durations_days[i] * 86400
        sigma = tau / (2 * np.sqrt(2 * np.log(2)))

        # perfil gaussiano
        burst_flux = (E_acc / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-(time - t0)**2 / (2 * sigma**2))
        lc += burst_flux

        burst_data.append({
            't0_day': t0 / 86400,
            'm_clump_Msun': m_clump / Msun,
            'E_acc_erg': E_acc,
            'tau_day': tau / 86400
        })

    # Extraer duraciones en días de los bursts
    durations_days = [b['tau_day'] for b in burst_data]

    # Definir nuevamente time_days para evitar el error
    time_days = time / 86400

    # Calcular tasa de acreción (M_dot)
    Mdot = 2 * lc / v_ff**2  # g/s
    Mdot_Msun_per_yr = Mdot / Msun * (3600 * 24 * 365.25)  # M_sun/yr

    # Añadir componente base de acreción y luminosidad correspondiente
    Mdot_Msun_per_yr += 1e-9
    # Calcular luminosidad base correspondiente
    Mdot_base = 1e-9 * Msun / (365.25 * 24 * 3600)  # g/s
    L_base = 0.5 * Mdot_base * v_ff**2  # erg/s
    lc += L_base

    t0_days = [b['t0_day'] for b in burst_data]
    masses_Msun = [b['m_clump_Msun'] for b in burst_data]

    return {
        'lc_model': lc,
        'time_days': time_days,
        'durations': durations_days,
        'masses': masses_Msun,
        't0_days': t0_days
    }

if __name__ == "__main__":
    results = run_clumpy_model()
    # Graficar curva de luz y tasa de acreción
    lc = results['lc_model']
    time_days = results['time_days']
    durations_days = results['durations']
    masses_Msun = results['masses']
    t0_days = results['t0_days']

    # Calcular tasa de acreción (M_dot) para graficar
    v_ff = np.sqrt(2 * G * M_star * (1/R_star - 1/R_trunc))
    Mdot = 2 * lc / v_ff**2  # g/s
    Mdot_Msun_per_yr = Mdot / Msun * (3600 * 24 * 365.25)  # M_sun/yr
    Mdot_Msun_per_yr += 1e-9

    fig, ax1 = plt.subplots(figsize=(12, 5))

    color = 'tab:blue'
    ax1.set_xlabel('Time (days)')
    #ax1.set_ylabel('Luminosity (10^33 erg/s)', color=color)
    #ax1.plot(time_days, lc / 1e33, color=color)
    #ax1.tick_params(axis='y', labelcolor=color)
    #ax1.grid(True)

    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Accretion Rate (M_sun/yr)', color=color)
    ax2.plot(time_days, Mdot_Msun_per_yr, color=color, alpha=0.6)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_yscale('log')

    plt.title('Clumpy Accretion: Luminosity and Accretion Rate for TW Hya (Modelo de Juguete)')
    plt.tight_layout()
    plt.show()

    # Graficar histogramas
    fig, axs = plt.subplots(1, 4, figsize=(21, 5))

    axs[0].hist(t0_days, bins=15, color='skyblue', edgecolor='black')
    axs[0].set_xlabel('Burst Time (days)')
    axs[0].set_ylabel('Number of Bursts')
    axs[0].set_title('Distribution of Burst Times')

    axs[1].hist(masses_Msun, bins=np.logspace(np.log10(min(masses_Msun)), np.log10(max(masses_Msun)), 12),
                density=True, color='salmon', edgecolor='black')
    axs[1].set_xscale('log')
    axs[1].set_xlabel(r'Clump Mass ($M_\odot$)')
    axs[1].set_ylabel(r'Density (1/$M_\odot$)')
    axs[1].set_title('Distribution of Clump Masses')

    # Histograma de duración
    axs[2].hist(durations_days, bins=np.logspace(np.log10(min(durations_days)), np.log10(max(durations_days)), 12),
                density=True, color='mediumseagreen', edgecolor='black')
    axs[2].set_xscale('log')
    axs[2].set_xlabel('Burst Duration (days)')
    axs[2].set_ylabel('Density (1/d)')
    axs[2].set_title('Distribution of Burst Durations')

    # Diagrama de dispersión: Clump Mass vs Burst Duration
    axs[3].scatter(masses_Msun, durations_days, color='darkorange', edgecolors='black', alpha=0.7)
    axs[3].set_xscale('log')
    axs[3].set_yscale('log')
    axs[3].set_xlabel(r'Clump Mass ($M_\odot$)')
    axs[3].set_ylabel('Burst Duration (days)')
    axs[3].set_title('Clump Mass vs Burst Duration')
    axs[3].grid(True, which="both", ls="--", linewidth=0.5)

    plt.tight_layout()
    plt.show()