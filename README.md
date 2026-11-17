# Paléo-acoustique — Modélisation Aéroacoustique du *Parasaurolophus*

Projet de recherche en **paléo-acoustique computationnelle** visant à reconstruire la signature acoustique et les vocalisations probables du *Parasaurolophus* (spécimen juvénile **RAM 14000**) à partir de modèles 3D fossiles (STL), d'extraction géométrique de conduits aériens et de synthèses biomécaniques.

---

## Architecture du Pipeline

### Solveurs scientifiques (existants)

```text
1. Acquisition 3D (Sketchfab API)          ──► scrapping.py
2. Extraction géométrique (Voxel + 3D Skel) ──► dataset.py / data-analysis.ipynb
   └── Génère : centerline_coordinates.npy & airway_features.csv
3. Analyse d'impédance acoustique (TMM)     ──► acoustic_reasonnance.py
4. Synthèse vocale & biomécanique           ──► 2_masses_larynx.py, rigorouscall.py,
                                                longrigorouscall.py, convolution_organique.py,
                                                gerative_odular.py, modularity_of_the_sound.py
```

### Orchestration multi-agents (Grokbot Paleo × Antigravity)

```text
État 1 HYPOTHESIS
  Grokbot lit models/<espece>/metadata.json
  → prompt d'ingénierie pour Antigravity (output/<run_id>/antigravity_prompt.md)

État 2 BUILD_AND_SIM
  Antigravity écrit/exécute le Python (trimesh / TMM / glotte ; Gmsh si besoin)
  → output/<run_id>/{source.wav, filtered.wav, metrics.json, spectrogram.png}

État 3 PEER_REVIEW
  Grokbot « Relecteur Critique » audite metrics + spectrogramme
  → cohérence f=c/4L ou c/2L, clipping, artefacts numériques

RETRY_LOOP
  Si échec → feedback à Antigravity (max 3 itérations)
```

Taxon seed : `models/parasaurolophus_ram14000/metadata.json`

```powershell
# Tests unitaires (M1–M7, < 15 s)
python -m pytest -q tests

# Mesures audio indépendantes (M1)
python -m pipeline.measure --run-dir output\parasaurolophus_ram14000_20261115T105317Z_eacbdf

# Solveur TMM Kirchhoff avec pertes (M2)
python -m pipeline.tmm_lossy --profile airway_features.csv --fmax 4000 --out output\tmm_check

# Synthèse couplée Kelly–Lochbaum + 2 masses (M4) et review non-circulaire (M3)
python -m pipeline.build_and_sim --taxon parasaurolophus_ram14000 --synth kl
python -m pipeline.build_and_sim --taxon parasaurolophus_ram14000 --synth legacy_iir

# Smoke orchestrateur local / mode agents
python -m pipeline.orchestrator --taxon parasaurolophus_ram14000 --mode local
python -m pipeline.orchestrator --taxon parasaurolophus_ram14000 --mode agents
python -m pipeline.orchestrator --taxon parasaurolophus_ram14000 --run-id <RUN_ID> --resume-review
```

Briefs agents : `prompts/grokbot_hypothesis.md`, `prompts/antigravity_build_sim.md`, `prompts/grokbot_peer_review.md`, `prompts/t3_agent_cycle.md`.

---

## Contenu du Dépôt

### Scripts Python & Notebooks

| Fichier | Rôle |
| :--- | :--- |
| `scrapping.py` | Recherche et téléchargement automatique de modèles 3D (*Parasaurolophus skull*) via l'API Sketchfab, décompression récursive des archives `.zip` imbriquées et extraction des `.stl`. |
| `dataset.py` | Chargement du maillage STL, voxelisation, squelettisation 3D (`skimage.morphology.skeletonize`), génération d'une ligne médiane ordonnée (`centerline_coordinates.npy`) et calcul des sections transversales (`area_mm2`, `perimeter_mm`, `hydraulic_diameter_mm`, `distance_mm`) exportées dans `airway_features.csv`. |
| `data-analysis.ipynb` | Notebook Jupyter d'analyse exploratoire : validation théorique (POC sur cylindre parfait), extraction sur le fossile et tracé des profils d'aire de section $A(s)$ et de diamètre hydraulique $D_h(s)$. |
| `acoustic_reasonnance.py` | Calcul du spectre d'impédance acoustique par **Méthode des Matrices de Transfert (TMM)** le long du conduit et identification des fréquences de résonance (formants $F_1, F_2, F_3$). |
| `2_masses_larynx.py` | Modèle biomécanique glottique à **2 masses couplées** (Steinecke & Herzel) résolu par Runge-Kutta (`RK45`) pour simuler la source laryngée basse fréquence. |
| `rigorouscall.py` | Simulation glottique stochastique (jitter neuromusculaire, collision tissulaire, bruit de turbulence) filtrée par les formants issus de l'analyse TMM. |
| `longrigorouscall.py` | Extension du modèle stochastique sur un cycle respiratoire complet (8 s) avec modulation musculaire basse fréquence. |
| `convolution_organique.py` | Génération de la **Réponse Impulsionnelle (IR)** du crâne fossile à partir des formants et convolution avec une source audio organique. |
| `gerative_odular.py` | Synthèse soustractive simple (onde en dent de scie + filtre IIR *peak* sur le premier formant + enveloppe d'attaque/relâchement). |
| `modularity_of_the_sound.py` | Synthèse multi-formantique dynamique avec glissando, vibrato et bruit de souffle filtré. |

### Données & Documentation

- `paleo_acoustic_data/RightSkeleton_RAM14000_printable2.stl` (~19 Mo) : Modèle 3D de base utilisé pour l'extraction géométrique.
- `centerline_coordinates.npy` : Coordonnées 3D de la ligne médiane (100 nœuds le long de l'axe antéro-postérieur).
- `airway_features.csv` : Tableau des propriétés géométriques et aéroacoustiques sur 99 tranches transversales.
- `ailog.md` : Journal technique des corrections et validations effectuées sur le pipeline d'extraction.
- `paleo_acoustic_data/Figure_1.png`, `Figure_2.png` : Figures générées lors de l'analyse.

*(Note : les fichiers audio `.wav` générés ne sont pas versionnés dans ce dépôt.)*

---

## Prérequis & Installation

```bash
pip install numpy pandas scipy trimesh scikit-image matplotlib librosa soundfile requests
```

### Utilisation rapide

1. **Extraction géométrique** (génère `centerline_coordinates.npy` et `airway_features.csv`) :
   ```bash
   python dataset.py
   ```

2. **Calcul des formants par Matrice de Transfert (TMM)** :
   ```bash
   python acoustic_reasonnance.py
   ```

3. **Synthèse d'un appel** :
   ```bash
   python rigorouscall.py
   ```

4. *(Optionnel)* **Téléchargement de spécimens Sketchfab** :
   Définir votre jeton API dans la variable d'environnement `SKETCHFAB_API_TOKEN` avant d'exécuter `scrapping.py` :
   ```powershell
   $env:SKETCHFAB_API_TOKEN = "votre_token_ici"
   python scrapping.py
   ```

