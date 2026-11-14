import numpy as np
from scipy.signal import iirpeak, lfilter, sawtooth, butter
from scipy.io import wavfile

def synthesize_advanced_dino_call(formants=[359.6, 1080.0, 2150.0], duration=4.0, sample_rate=44100):
    print("Génération du cri organique avec dynamique non-linéaire...")
    t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    
    # ---------------------------------------------------------
    # STRATÉGIE 1 : Modulation Dynamique (Glissando + Vibrato)
    # ---------------------------------------------------------
    # Profil de l'appel : commence à 45Hz, monte à 75Hz, redescend à 50Hz
    f0_contour = np.interp(t, [0, duration*0.3, duration], [45.0, 75.0, 50.0])
    
    # Ajout d'un vibrato organique (instabilité de 6 Hz avec une amplitude de 2 Hz)
    vibrato = 2.0 * np.sin(2 * np.pi * 6.0 * t)
    f0_dynamic = f0_contour + vibrato
    
    # Calcul de la phase instantanée (intégrale de la fréquence)
    phase = 2 * np.pi * np.cumsum(f0_dynamic) / sample_rate
    
    # Source harmonique
    harmonic_source = sawtooth(phase)
    
    # ---------------------------------------------------------
    # STRATÉGIE 2 : Turbulence (Bruit de souffle)
    # ---------------------------------------------------------
    # Génération d'un bruit blanc
    noise = np.random.normal(0, 0.5, len(t))
    
    # Filtre passe-bas sur le bruit pour simuler les tissus charnus (coupure à 800 Hz)
    b_lp, a_lp = butter(2, 800.0 / (sample_rate / 2), btype='low')
    breath_noise = lfilter(b_lp, a_lp, noise)
    
    # Enveloppe dynamique du bruit (plus de souffle au début de l'appel)
    noise_envelope = np.interp(t, [0, 0.5, duration], [1.0, 0.3, 0.1])
    
    # Mixage : 85% cordes vocales, 15% turbulence aéroacoustique variable
    mixed_source = (harmonic_source * 0.85) + (breath_noise * noise_envelope * 0.15)
    
    # ---------------------------------------------------------
    # STRATÉGIE 3 : Filtrage Multi-Formantique
    # ---------------------------------------------------------
    filtered_signal = mixed_source.copy()
    
    # On applique les résonances du crâne en cascade
    # Les formants supérieurs ont généralement un facteur Q plus élevé (pics plus étroits) et moins d'amplitude
    gains = [1.0, 0.6, 0.3] 
    
    final_resonances = np.zeros_like(t)
    for f, gain in zip(formants, gains):
        b_peak, a_peak = iirpeak(f, Q=10.0, fs=sample_rate)
        final_resonances += lfilter(b_peak, a_peak, filtered_signal) * gain

    # ---------------------------------------------------------
    # Enveloppe finale et Export
    # ---------------------------------------------------------
    attack_samples = int(0.6 * sample_rate)
    release_samples = int(1.5 * sample_rate)
    
    envelope = np.ones_like(t)
    envelope[:attack_samples] = np.linspace(0, 1, attack_samples)**2 # Courbe d'attaque exponentielle
    envelope[-release_samples:] = np.linspace(1, 0, release_samples)**2
    
    final_audio = final_resonances * envelope
    
    # Normalisation
    final_audio = final_audio / np.max(np.abs(final_audio))
    audio_int16 = np.int16(final_audio * 32767)
    
    wavfile.write("parasaurolophus_organic_call.wav", sample_rate, audio_int16)
    print("Fichier 'parasaurolophus_organic_call.wav' généré avec succès.")

if __name__ == '__main__':
    # Remplacer formants par les F1, F2, F3 complets trouvés par l'algorithme TMM
    synthesize_advanced_dino_call(formants=[359.6, 1080.0, 2150.0])