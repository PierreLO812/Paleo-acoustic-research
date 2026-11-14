import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def compute_acoustic_resonance(csv_path, max_freq=1000, num_points=2000):
    # 1. Chargement des données géométriques
    df = pd.read_csv(csv_path)
    
    # Conversion en unités SI (Mètres et Mètres carrés)
    areas = df['area_mm2'].values * 1e-6
    distances = df['distance_mm'].values * 1e-3
    
    # Calcul de la longueur (L) et de l'aire moyenne (S) de chaque segment (tranche)
    L = np.diff(distances)
    S = (areas[:-1] + areas[1:]) / 2
    
    # 2. Constantes physiques (Air à ~20°C, modifiable pour un dinosaure exothermique/endothermique)
    c = 343.0      # Vitesse du son (m/s)
    rho = 1.225    # Densité de l'air (kg/m3)
    
    freqs = np.linspace(10, max_freq, num_points)
    Z_in = np.zeros(len(freqs), dtype=complex)
    
    print("Calcul de l'impédance acoustique (Matrice de Transfert)...")
    for i, f in enumerate(freqs):
        k = (2 * np.pi * f) / c
        
        # Condition aux limites aux narines : Z_L approx 0 (tube ouvert)
        # (Pour plus de précision, on pourrait ajouter l'impédance de radiation)
        Z = 0 + 0j 
        
        # Remontée du conduit des narines vers la glotte
        for j in range(len(L)-1, -1, -1):
            Z_c = rho * c / S[j] # Impédance caractéristique du segment
            
            # Équation de transfert d'impédance
            num = Z + 1j * Z_c * np.tan(k * L[j])
            den = Z_c + 1j * Z * np.tan(k * L[j])
            Z = Z_c * (num / den)
            
        Z_in[i] = Z

    # 3. Identification des pics de résonance (Formants)
    Z_mag = np.abs(Z_in)
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(Z_mag)
    resonance_freqs = freqs[peaks]
    
    print(f"Fréquences de résonance fondamentales trouvées (F1, F2, F3) :")
    for i, rf in enumerate(resonance_freqs[:3]):
        print(f" - Formant {i+1} : {rf:.1f} Hz")
        
    # 4. Visualisation
    plt.figure(figsize=(10, 5))
    plt.plot(freqs, Z_mag, color='indigo', lw=2)
    plt.plot(freqs[peaks], Z_mag[peaks], "rx", label="Pics de résonance")
    plt.title("Spectre d'Impédance Acoustique du Parasaurolophus")
    plt.xlabel("Fréquence (Hz)")
    plt.ylabel("Magnitude de l'Impédance (Ohms acoustiques)")
    plt.xlim(0, max_freq)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    compute_acoustic_resonance('airway_features.csv')