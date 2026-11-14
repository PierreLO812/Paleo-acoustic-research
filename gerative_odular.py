import numpy as np
from scipy.signal import iirpeak, lfilter, sawtooth
from scipy.io import wavfile

def synthesize_dino_call(formant_freq, f0=60.0, duration=4.0, sample_rate=44100):
    """
    Génère un cri de dinosaure basé sur la fréquence de résonance du crâne.
    f0 : Fréquence fondamentale (cordes vocales/syrinx). 60 Hz donne un grondement profond.
    formant_freq : La fréquence de résonance calculée (359.6 Hz).
    """
    print(f"Génération du signal audio (F0: {f0} Hz, Formant: {formant_freq} Hz)...")
    
    # 1. Axe du temps
    t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    
    # 2. Le Signal Source (L'organe vocal)
    # L'onde en dent de scie est idéale car elle contient naturellement tous les harmoniques pairs et impairs
    source = sawtooth(2 * np.pi * f0 * t)
    
    # 3. Le Filtre (Le crâne fossile)
    # Création d'un filtre IIR "Peak" pour amplifier spécifiquement la résonance
    Q = 8.0  # Facteur de qualité (étroitesse du pic de résonance)
    b, a = iirpeak(formant_freq, Q, sample_rate)
    
    # Passage du son des cordes vocales à travers le crâne
    filtered_signal = lfilter(b, a, source)
    
    # 4. L'Enveloppe (Le souffle organique)
    # Un son constant fait "robotique". On applique une enveloppe d'amplitude (Fade-in / Hold / Fade-out)
    attack_time = 0.8  # secondes pour atteindre le volume max (inspiration/poussée)
    release_time = 1.5 # secondes pour que le son s'éteigne
    
    attack_samples = int(attack_time * sample_rate)
    release_samples = int(release_time * sample_rate)
    
    envelope = np.ones_like(t)
    envelope[:attack_samples] = np.linspace(0, 1, attack_samples)
    envelope[-release_samples:] = np.linspace(1, 0, release_samples)
    
    # Application de l'enveloppe
    final_audio = filtered_signal * envelope
    
    # 5. Normalisation et Export
    # Empêche la saturation (clipping) et convertit en PCM 16-bit standard
    final_audio = final_audio / np.max(np.abs(final_audio))
    audio_int16 = np.int16(final_audio * 32767)
    
    output_filename = "parasaurolophus_call.wav"
    wavfile.write(output_filename, sample_rate, audio_int16)
    print(f"Succès ! Fichier exporté : {output_filename}")

if __name__ == '__main__':
    # Injection du Formant 1 calculé à l'étape précédente
    synthesize_dino_call(formant_freq=359.6)