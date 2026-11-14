import numpy as np
import librosa
import soundfile as sf
from scipy.signal import convolve, get_window

def create_fossil_impulse_response(formants, bandwidths, duration=0.5, sr=44100):
    """
    Génère la Réponse Impulsionnelle (IR) du crâne fossile.
    C'est la signature acoustique exacte de la cavité 3D.
    """
    t = np.linspace(0, duration, int(sr * duration), endpoint=False)
    ir = np.zeros_like(t)
    
    # Pour chaque formant de ton analyse TMM, on génère une résonance qui s'estompe (damped sine)
    for f, bw in zip(formants, bandwidths):
        # Le taux de déclin (decay) est lié à la bande passante (l'étroitesse du pic)
        decay = np.exp(-np.pi * bw * t)
        resonance = np.sin(2 * np.pi * f * t) * decay
        ir += resonance
        
    # On applique une fenêtre pour lisser la fin de la réverbération interne
    ir *= get_window('hann', len(ir))
    # Normalisation de l'IR
    return ir / np.max(np.abs(ir))

def convolve_fossil_audio(organic_audio_path, output_path, formants, bandwidths):
    print(f"Chargement de la source organique : {organic_audio_path}")
    
    # 1. Charger le son animal réel (ex: alligator.wav)
    # librosa convertit automatiquement en mono et gère les fréquences d'échantillonnage
    source_audio, sr = librosa.load(organic_audio_path, sr=44100)
    
    # 2. Générer la géométrie acoustique du fossile
    print("Génération de l'Impulse Response du Parasaurolophus...")
    fossil_ir = create_fossil_impulse_response(formants, bandwidths, sr=sr)
    
    # 3. CONVOLUTION (La magie opère ici)
    # Le son de l'animal est "projeté" à l'intérieur du crâne fossile
    print("Application de la convolution aéroacoustique...")
    convolved_audio = convolve(source_audio, fossil_ir, mode='full')
    
    # 4. Nettoyage et Export
    # Un léger filtre passe-haut pour enlever les infra-basses inutiles créées par la convolution
    convolved_audio = librosa.effects.preemphasis(convolved_audio, coef=0.85)
    
    # Normalisation stricte pour éviter la saturation
    convolved_audio = convolved_audio / np.max(np.abs(convolved_audio))
    
    sf.write(output_path, convolved_audio, sr)
    print(f"Rendu biomécanique sauvegardé : {output_path}")

if __name__ == '__main__':
    # Fréquences issues de ton analyse (F1, F2, F3)
    # Les "bandwidths" définissent la largeur des pics d'impédance (en Hz). Plus c'est petit, plus ça résonne longtemps.
    fossil_formants = [359.6, 1080.0, 2150.0]
    fossil_bandwidths = [40.0, 90.0, 150.0] 
    
    # Il te faut un fichier source réel pour tester (ex: un alligator ou un éléphant)
convolve_fossil_audio('biomechanical_larynx.wav', 'parasaurolophus_organic_call.wav', fossil_formants, fossil_bandwidths)