# Antigravity — BUILD_AND_SIM for `parasaurolophus_ram14000`

You are the engineering agent for Pierre's Paleo-acoustic research repo.
Repo root: `C:\Users\balad\OneDrive\Documents\Paléo-acoustique\paleoacoustique`
Run id: `parasaurolophus_ram14000_20261117T223414Z_2b13b0`
Output directory (create it): `C:\Users\balad\OneDrive\Documents\Paléo-acoustique\paleoacoustique\output\parasaurolophus_ram14000_20261117T223414Z_2b13b0`

## Scientific goal
Produce a physically constrained vocalization sketch for **Parasaurolophus sp.** (RAM14000),
not a cinematic roar. Soft tissue is hypothesized; geometry comes from fossil STL / extracted airway.

## Raw anatomical parameters (from models/parasaurolophus_ram14000/metadata.json)
```json
{
  "taxon": "Parasaurolophus sp.",
  "specimen_id": "RAM14000",
  "common_label": "parasaurolophus_ram14000",
  "clade": "Hadrosauridae",
  "life_stage": "juvenile",
  "geometry": {
    "stl_path": "paleo_acoustic_data/RightSkeleton_RAM14000_printable2.stl",
    "airway_features_csv": "airway_features.csv",
    "centerline_npy": "centerline_coordinates.npy",
    "conduit_length_mm": 263.804,
    "mean_area_mm2": 1371.2,
    "mean_hydraulic_diameter_mm": 20.49,
    "notes": "Centerline and cross-sections from trimesh voxelization + skimage skeletonize (dataset.py)."
  },
  "soft_tissue_assumptions": {
    "lining_present": true,
    "lining_thickness_mm": 2.0,
    "superseded_by": "soft_tissue_hypotheses.json",
    "glottal_source": "two_mass_stochastic",
    "boundary_condition": "open_nostril_approx_Z0",
    "temperature_c": 20.0,
    "speed_of_sound_m_s": 343.0,
    "confidence": 0.55,
    "notes": "Soft tissues not preserved; lining thickness and glottal parameters are explicit priors."
  },
  "acoustic_priors": {
    "tube_model": "closed_glottis_open_nostril",
    "theoretical_formula": "f_n = (2n-1) * c / (4L)",
    "expected_formants_hz": [
      359.6,
      1080.0,
      2150.0
    ],
    "default_q_factors": [
      5.0,
      8.0,
      10.0
    ],
    "default_bandwidths_hz": [
      70.0,
      135.0,
      215.0
    ],
    "f0_prior_hz": 60.0
  },
  "solvers": {
    "geometry": "dataset.py",
    "tmm_impedance": "acoustic_reasonnance.py",
    "glottal_source": "rigorouscall.py",
    "optional_long_call": "longrigorouscall.py",
    "mesh_stack": [
      "trimesh",
      "skimage.morphology.skeletonize"
    ],
    "future_mesh": [
      "gmsh"
    ]
  },
  "pipeline": {
    "required_artifacts": [
      "source.wav",
      "filtered.wav",
      "metrics.json",
      "spectrogram.png"
    ],
    "max_retries": 3
  }
}
```

## Theoretical resonance anchors
- Conduit length L = 0.263804 m
- c = 343.0 m/s
- Closed-open quarters: [325.05, 975.16, 1625.26] Hz  (f = (2n-1)c / 4L)
- Open-open halves: [650.1, 1300.21, 1950.31] Hz  (f = n c / 2L)
- Preferred model: `f_n = (2n-1) * c / (4L)`
- Prior formants from metadata: [359.6, 1080.0, 2150.0] Hz, Q≈[5.0, 8.0, 10.0] (targets only, never copy as measured peaks)

## Mandatory engineering steps
1. Use the independent solvers:
   - Geometry QC: `pipeline/geometry_qc.py`
   - Lossy Kirchhoff TMM: `pipeline/tmm_lossy.py`
   - Coupled KL waveguide + two-mass source: `pipeline/source_tract.py`
   - Audio-only measurement: `pipeline/measure.py`
2. Call `python -m pipeline.build_and_sim --taxon parasaurolophus_ram14000 --run-id parasaurolophus_ram14000_20261117T223414Z_2b13b0`.
3. STL reference: `paleo_acoustic_data/RightSkeleton_RAM14000_printable2.stl` ; airway table: `airway_features.csv`
4. Do **not** invent fossil measurements. If a file is missing, fail clearly.

## Mandatory artifacts in `C:\Users\balad\OneDrive\Documents\Paléo-acoustique\paleoacoustique\output\parasaurolophus_ram14000_20261117T223414Z_2b13b0/`
- `source.wav`, `filtered.wav` (-1 dBFS peak normalized, no clipping)
- `spectrogram.png`, `transfer_function.png`, `tmm_lossy.png`, `kl_vs_tmm.png`, `geometry_qc.png`
- `metrics.json` (schema v2.0), `measured.json`, `tmm_lossy.json`, `state_*.json`, `report.md`, `build_log.txt`


## Done criteria
All artifacts exist, `metrics.json` validates against `schemas/metrics.schema.json`, and `report.md` / `build_log.txt` are written.
