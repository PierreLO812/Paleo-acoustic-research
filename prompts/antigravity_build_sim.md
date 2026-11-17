# Antigravity — État BUILD_AND_SIM

Tu reçois un prompt généré en HYPOTHESIS (`output/<run_id>/antigravity_prompt.md`).

## Mission
1. Écrire / adapter le Python en t'appuyant sur les solveurs existants (trimesh/`dataset.py`, TMM/`acoustic_reasonnance.py`, glotte/`rigorouscall.py`; Gmsh si nécessaire).
2. Exécuter le calcul.
3. Produire **exactement** dans `output/<run_id>/` :
   - `source.wav`
   - `filtered.wav`
   - `metrics.json` (pics de formants Hz, bandes passantes, facteurs Q, bloc théorique, clipping)
   - `spectrogram.png`

## Chemin rapide valide
```powershell
python -m pipeline.build_and_sim --taxon parasaurolophus_ram14000 --run-id <RUN_ID>
```
Améliore ce module si le peer-review échoue (ne casse pas le contrat d'artefacts).

## Interdits
- Inventer des mesures fossiles absentes de `metadata.json`
- Remplacer la physique par un effet "cinéma"
- Omettre un artefact obligatoire
