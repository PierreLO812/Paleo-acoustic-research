import numpy as np
from scipy.integrate import solve_ivp
from scipy.io import wavfile
import matplotlib.pyplot as plt

def vocal_fold_biomechanics(t, state, params):
    """
    Système d'équations différentielles du modèle à 2 masses (Steinecke & Herzel).
    state = [x1, v1, x2, v2] (positions et vélocités des masses 1 et 2)
    """
    x1, v1, x2, v2 = state
    m1, m2, k1, k2, kc, r1, r2, P_sub, c_col = params
    
    # 1. Aire glottique (ouverture entre les cordes)
    # L'aire au repos est a0. Si x < -a0, il y a collision.
    a0 = 0.05 
    a1 = a0 + x1
    a2 = a0 + x2
    
    # 2. Forces de collision (non-linéarité extrême simulant le claquement des tissus)
    # Si l'aire devient négative, les tissus s'écrasent, créant une force de répulsion massive
    F_col1 = -c_col * a1 if a1 < 0 else 0
    F_col2 = -c_col * a2 if a2 < 0 else 0
    
    # 3. Forces aérodynamiques (Pression sous-glottique et effet Bernoulli)
    # Modèle simplifié de la pression sur les masses
    # La pression s'applique surtout sur m1 (masse inférieure) quand la glotte est ouverte
    if a1 > 0 and a2 > 0:
        # L'air passe : Bernoulli accélère le flux et crée une dépression (qui referme les cordes)
        F_aero1 = P_sub * (a0 / a1)**2
        F_aero2 = 0.0  # Simplification : dissipation de la pression sur la masse supérieure
    elif a1 > 0 and a2 <= 0:
        # m1 ouverte, m2 fermée : P_sub s'applique à fond sur m1 pour forcer l'ouverture
        F_aero1 = P_sub
        F_aero2 = 0.0
    else:
        # Totalement fermé
        F_aero1 = P_sub
        F_aero2 = 0.0

    # 4. Équations du mouvement (Oscillateurs harmoniques amortis, couplés et forcés)
    # Somme des forces = - Rappel élastique - Amortissement - Couplage + Aérodynamique + Collision
    dx1_dt = v1
    dv1_dt = (1/m1) * (-k1*x1 - r1*v1 - kc*(x1 - x2) + F_aero1 + F_col1)
    
    dx2_dt = v2
    dv2_dt = (1/m2) * (-k2*x2 - r2*v2 - kc*(x2 - x1) + F_aero2 + F_col2)
    
    return [dx1_dt, dv1_dt, dx2_dt, dv2_dt]

def generate_biomechanical_source(duration=1.0, sr=44100):
    print("Résolution des équations biomécaniques (RK45)...")
    
    # PARAMÈTRES PHYSIQUES (Ajustés pour un grand animal type dinosaure)
    # Les masses sont plus lourdes que l'humain pour abaisser la fréquence (F0 ~ 40-60 Hz)
    m1, m2 = 0.5, 0.3        # Masses (g)
    k1, k2 = 40.0, 30.0      # Raideur des ressorts longitudinaux (tension)
    kc = 15.0                # Couplage entre les deux masses (rigidité du tissu)
    r1, r2 = 1.2, 0.8        # Amortissement (viscosité des tissus)
    P_sub = 80.0             # Pression sous-glottique (la puissance du souffle)
    c_col = 300.0            # Raideur de collision (plus c'est élevé, plus le son est "rocailleux")
    
    params = (m1, m2, k1, k2, kc, r1, r2, P_sub, c_col)
    
    # Conditions initiales (léger décalage pour amorcer le chaos)
    y0 = [0.001, 0.0, 0.0, 0.0] 
    
    # Vecteur temps précis pour correspondre à la fréquence d'échantillonnage audio
    t_span = (0, duration)
    t_eval = np.linspace(0, duration, int(sr * duration))
    
    # Résolution du système via Runge-Kutta
    sol = solve_ivp(vocal_fold_biomechanics, t_span, y0, t_eval=t_eval, args=(params,), method='RK45')
    
    # Le flux d'air acoustique (Glottal Volume Velocity) est proportionnel à l'ouverture minimale
    # On prend la dérivée de l'aire pour simuler la pression acoustique rayonnée
    glottal_area = np.maximum(0, np.minimum(sol.y[0] + 0.05, sol.y[2] + 0.05))
    glottal_flow = np.gradient(glottal_area)
    
    # Normalisation pour l'audio
    audio_signal = glottal_flow / np.max(np.abs(glottal_flow))
    
    return t_eval, glottal_area, audio_signal, sr

if __name__ == '__main__':
    # 1. Génération du signal
    t, area, audio, sample_rate = generate_biomechanical_source(duration=2.0)
    
    # 2. Export du fichier audio source (L'excitateur)
    wavfile.write("biomechanical_larynx.wav", sample_rate, np.int16(audio * 32767))
    print("Source glottique pure exportée : 'biomechanical_larynx.wav'")
    
    # 3. Visualisation de la vibration chaotique (sur 100 ms)
    mask = t < 0.1
    plt.figure(figsize=(10, 4))
    plt.plot(t[mask], area[mask], color='darkred', lw=1.5)
    plt.title("Dynamique chaotique de l'aire glottique (Modèle à 2 masses)")
    plt.xlabel("Temps (s)")
    plt.ylabel("Ouverture glottique")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()