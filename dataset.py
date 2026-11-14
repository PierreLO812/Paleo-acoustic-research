import os
import sys
import numpy as np
import pandas as pd
import trimesh

# Compatibilité scikit-image (skeletonize_3d a été remplacé par skeletonize)
try:
    from skimage.morphology import skeletonize
except ImportError:
    from skimage.morphology import skeletonize_3d as skeletonize

def find_stl_file(candidate_name="RightSkeleton_RAM14000_printable2.stl"):
    """Localise le fichier STL dans le répertoire courant ou dans paleo_acoustic_data."""
    # 1. Vérification par argument en ligne de commande
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        return sys.argv[1]

    # 2. Chemins probables standards
    direct_candidates = [
        candidate_name,
        os.path.join("paleo_acoustic_data", candidate_name),
        os.path.join(".", "paleo_acoustic_data", candidate_name),
    ]
    for c in direct_candidates:
        if os.path.exists(c):
            return os.path.abspath(c)

    # 3. Recherche récursive dans paleo_acoustic_data et le répertoire courant
    search_dirs = ["paleo_acoustic_data", "."]
    for sdir in search_dirs:
        if os.path.exists(sdir):
            for root, _, files in os.walk(sdir):
                for f in files:
                    if f.lower() == candidate_name.lower():
                        return os.path.abspath(os.path.join(root, f))

    # 4. Fallback vers n'importe quel fichier STL trouvé
    for sdir in search_dirs:
        if os.path.exists(sdir):
            for root, _, files in os.walk(sdir):
                for f in files:
                    if f.lower().endswith(".stl"):
                        return os.path.abspath(os.path.join(root, f))

    raise FileNotFoundError(f"Fichier STL introuvable : '{candidate_name}'. Vérifiez son emplacement.")

def main():
    stl_path = find_stl_file("RightSkeleton_RAM14000_printable2.stl")
    print(f"Fichier STL sélectionné : {stl_path}")

    # 1. CHARGEMENT DU MAILLAGE
    print("1. Chargement du STL...")
    mesh = trimesh.load(stl_path, force='mesh')
    print(f"   Maillage chargé : {len(mesh.vertices)} sommets, {len(mesh.faces)} faces.")
    print(f"   Dimensions (extents) : {np.round(mesh.extents, 2)} mm")

    # 2. VOXELISATION
    print("2. Voxelisation en cours...")
    pitch = mesh.extents.max() / 100.0
    print(f"   Résolution de voxel (pitch) : {pitch:.3f} mm")
    voxels = mesh.voxelized(pitch=pitch)

    # 3. SQUELETTISATION (EXTRACTION DE LA LIGNE MÉDIANE)
    print("3. Extraction du squelette (skeletonization 3D)...")
    skeleton_matrix = skeletonize(voxels.matrix)
    indices = np.argwhere(skeleton_matrix)
    raw_pts = voxels.indices_to_points(indices)
    print(f"   {len(raw_pts)} voxels squelettes identifiés.")

    # 4. AGRÉGATION ET LISSAGE LE LONG DU CONDUIT
    axe_dominant = np.argmax(mesh.extents)
    axes_noms = ['X', 'Y', 'Z']
    print(f"   Axe dominant détecté : {axes_noms[axe_dominant]}")

    # Regroupement des voxels par tranche le long de l'axe dominant pour obtenir une vraie ligne médiane
    unique_coords, inverse_indices = np.unique(raw_pts[:, axe_dominant], return_inverse=True)
    centerline_points = np.array([raw_pts[inverse_indices == i].mean(axis=0) for i in range(len(unique_coords))])
    centerline_points = centerline_points[np.argsort(centerline_points[:, axe_dominant])]

    # Sauvegarde des coordonnées de la ligne médiane pour réutilisation directe
    np.save('centerline_coordinates.npy', centerline_points)
    print(f"4. Ligne médiane générée : {len(centerline_points)} nœuds (sauvegardée dans 'centerline_coordinates.npy').")

    # 5. CALCUL DES SECTIONS TRANSVERSALES ET DES PROPRIÉTÉS AÉROACOUSTIQUES
    print("5. Calcul des sections transversales et métriques aéroacoustiques...")
    dataset = []

    for i in range(len(centerline_points) - 1):
        p1 = centerline_points[i]
        p2 = centerline_points[i + 1]

        # Vecteur tangent le long du flux d'air
        direction = p2 - p1
        norm = np.linalg.norm(direction)
        if norm == 0:
            continue
        direction = direction / norm

        # Distance cumulative le long de la ligne médiane
        dist = np.linalg.norm(centerline_points[0] - p1) if i > 0 else 0.0

        # Coupe géométrique du STL
        try:
            slice_2d = mesh.section(plane_origin=p1, plane_normal=direction)
        except Exception:
            continue

        area, perimeter = 0.0, 0.0

        if slice_2d is not None:
            try:
                # Utilisation de to_2D() (remplace to_planar() déprécié dans Trimesh)
                if hasattr(slice_2d, 'to_2D'):
                    planar_slice, _ = slice_2d.to_2D()
                else:
                    planar_slice, _ = slice_2d.to_planar()
                area = planar_slice.area
                perimeter = planar_slice.length
            except Exception:
                pass

        # Diamètre hydraulique Dh = 4 * Aire / Périmètre
        hydraulic_diam = (4 * area / perimeter) if perimeter > 0 else 0.0

        dataset.append({
            'node_id': i,
            'distance_mm': round(dist, 3),
            'x': round(p1[0], 3),
            'y': round(p1[1], 3),
            'z': round(p1[2], 3),
            'area_mm2': round(area, 3),
            'perimeter_mm': round(perimeter, 3),
            'hydraulic_diameter_mm': round(hydraulic_diam, 3)
        })

    # 6. EXPORT EN DATASET
    df_airway = pd.DataFrame(dataset)
    df_airway.to_csv('airway_features.csv', index=False)
    print("6. Succès ! Dataset exporté sous 'airway_features.csv'")
    print(f"   Nombre total de tranches calculées : {len(df_airway)}")
    print(f"   Longueur totale approximée : {df_airway['distance_mm'].max():.1f} mm")
    print(f"   Aire moyenne de section : {df_airway['area_mm2'].mean():.2f} mm² (min: {df_airway['area_mm2'].min():.2f}, max: {df_airway['area_mm2'].max():.2f})")
    print(f"   Diamètre hydraulique moyen : {df_airway['hydraulic_diameter_mm'].mean():.2f} mm")

if __name__ == '__main__':
    main()