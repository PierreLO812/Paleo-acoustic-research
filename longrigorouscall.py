import numpy as np
import scipy.signal as signal
from scipy.io import wavfile

def generate_extended_dino_call(duration=8.0, sr=44100):
    print(f"1. Résolution de l'aérodynamique sur un cycle de {duration} secondes...")
    
    n_steps = int(duration * sr)
    dt = 1.0 / sr
    t = np.linspace(0, duration, n_steps)
    
    # Paramètres biomécaniques optimisés pour la collision
    m1, m2 = 0.9, 0.6
    k1_base, k2_base = 35.0, 18.0
    kc = 12.0
    r1, r2 = 0.3, 0.15
    c_col = 1000.0  # Choc tissulaire très dur pour un timbre agressif
    a0 = 0.002      # Glotte très resserrée
    
    # Enveloppe de pression pulmonaire (cycle respiratoire complexe)
    # 0s: début | 0.5s: attaque max | 1.5s: creux | 3.0s: reprise | 7.0s: essoufflement
    times = [0, 0.5, 1.5, 3.0, 5.0, 7.0, duration]
    pressures = [0, 260.0, 200.0, 240.0, 180.0, 80.0, 0]
    P_sub_env = np.interp(t, times, pressures)
    
    # Modulation musculaire basse fréquence (simule l'effort organique)
    muscle_mod = 1.0 + 0.08 * np.sin(2 * np.pi * 1.5 * t) + 0.04 * np.sin(2 * np.pi * 4.2 * t)
    
    # Pré-génération du Jitter (bruit neuromusculaire à haute fréquence)
    jitter1 = 1.0 + 0.03 * np.random.randn(n_steps)
    jitter2 = 1.0 + 0.03 * np.random.randn(n_steps)
    
    x1, v1 = np.zeros(n_steps), np.zeros(n_steps)
    x2, v2 = np.zeros(n_steps), np.zeros(n_steps)
    glottal_area = np.zeros(n_steps)
    
    x1[0], x2[0] = 0.01, -0.01 
    
    # Intégration d'Euler rapide
    for i in range(n_steps - 1):
        k1 = k1_base * muscle_mod[i] * jitter1[i]
        k2 = k2_base * muscle_mod[i] * jitter2[i]
        P = P_sub_env[i]
        
        a1 = a0 + x1[i]
        a2 = a0 + x2[i]
        
        F_col1 = -c_col * a1 if a1 < 0 else 0
        F_col2 = -c_col * a2 if a2 < 0 else 0
        
        if a1 > 0 and a2 > 0:
            F_aero1 = P * (a0 / a1)**2
        elif a1 > 0 and a2 <= 0:
            F_aero1 = P
        else:
            F_aero1 = P
            
        F_aero2 = 0.0
        
        a1_acc = (1/m1) * (-k1*x1[i] - r1*v1[i] - kc*(x1[i] - x2[i]) + F_aero1 + F_col1)
        a2_acc = (1/m2) * (-k2*x2[i] - r2*v2[i] - kc*(x2[i] - x1[i]) + F_aero2 + F_col2)
        
        v1[i+1] = v1[i] + a1_acc * dt
        x1[i+1] = x1[i] + v1[i+1] * dt
        
        v2[i+1] = v2[i] + a2_acc * dt
        x2[i+1] = x2[i] + v2[i+1] * dt
        
        glottal_area[i] = max(0, min(a0 + x1[i], a0 + x2[i]))
        
    glottal_flow = np.gradient(glottal_area)
    
    # Bruit d'aspiration minimal (limité à 1.5%)
    turbulence = np.random.normal(0, 1, n_steps)
    b_lp, a_lp = signal.butter(2, 800 / (sr / 2), btype='low')
    turbulence = signal.lfilter(b_lp, a_lp, turbulence)
    noise_mod = turbulence * glottal_area * (P_sub_env / np.max(P_sub_env)) * 0.015
    
    source_audio = glottal_flow + noise_mod

    print("2. Filtrage aéroacoustique par Matrice de Transfert (TMM)...")
    formants = [359.6, 1080.0, 2150.0]
    gains = [1.0, 0.45, 0.15]
    # Q_factors augmentés pour accentuer l'effet "caverneux" du crâne
    Q_factors = [6.0, 8.0, 12.0]
    
    final_audio = np.zeros_like(source_audio)
    for f, gain, q in zip(formants, gains, Q_factors):
        b, a = signal.iirpeak(f, Q=q, fs=sr)
        final_audio += signal.lfilter(b, a, source_audio) * gain
        
    # Lissage des extrémités pour éviter les "clics" numériques
    master_env = np.ones_like(t)
    fade_in = int(0.1 * sr)
    fade_out = int(0.5 * sr)
    master_env[:fade_in] = np.linspace(0, 1, fade_in)
    master_env[-fade_out:] = np.linspace(1, 0, fade_out)
    final_audio *= master_env
        
    final_audio = final_audio / np.max(np.abs(final_audio))
    output_path = "parasaurolophus_8s_call.wav"
    wavfile.write(output_path, sr, np.int16(final_audio * 32767))
    print(f"3. Terminé ! Fichier de {duration}s exporté sous : {output_path}")

if __name__ == '__main__':
    generate_extended_dino_call()