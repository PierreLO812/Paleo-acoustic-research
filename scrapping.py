import os
import re
import shutil
import zipfile
import requests

# Configuration (définir la variable d'environnement SKETCHFAB_API_TOKEN)
API_TOKEN = os.environ.get('SKETCHFAB_API_TOKEN', '')  # Token Sketchfab
SEARCH_QUERY = 'Parasaurolophus skull'
OUTPUT_DIR = os.path.abspath('./paleo_acoustic_data')

HEADERS = {
    'Authorization': f'Token {API_TOKEN}'
}

def sanitize_filename(name: str) -> str:
    """Remplace les caractères interdits dans les noms de fichiers Windows."""
    return re.sub(r'[\\/*?:"<>|]', '_', name).strip()

def extract_nested_zips(directory: str):
    """Parcourt récursivement le dossier pour décompresser les archives zip imbriquées."""
    found_nested = True
    while found_nested:
        found_nested = False
        for root, _, files in os.walk(directory):
            for file in files:
                if file.lower().endswith('.zip'):
                    nested_zip_path = os.path.join(root, file)
                    target_dir = os.path.splitext(nested_zip_path)[0]
                    # Évite d'extraire indéfiniment si le dossier cible existe déjà
                    if not os.path.exists(target_dir):
                        print(f"Extraction de l'archive imbriquée : {nested_zip_path}")
                        try:
                            with zipfile.ZipFile(nested_zip_path, 'r') as zf:
                                zf.extractall(target_dir)
                            found_nested = True
                        except zipfile.BadZipFile:
                            print(f"Archive corrompue : {nested_zip_path}")

def find_and_expose_stl_files(directory: str, output_directory: str):
    """Identifie les fichiers STL extraits et les rend facilement accessibles."""
    stl_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith('.stl'):
                full_path = os.path.join(root, file)
                stl_files.append(full_path)
                # Si le fichier n'est pas déjà dans le dossier racine de sortie, on en place une copie directe
                direct_path = os.path.join(output_directory, file)
                if not os.path.exists(direct_path) and os.path.abspath(full_path) != os.path.abspath(direct_path):
                    try:
                        shutil.copy2(full_path, direct_path)
                        print(f"Copie pratique vers le dossier racine : {direct_path}")
                    except Exception as e:
                        print(f"Impossible de copier {file} : {e}")

    print(f"\nFichiers STL disponibles ({len(stl_files)}) :")
    for stl in stl_files:
        print(f" - {stl}")
    return stl_files

def search_and_download_specimens(query: str, output_directory: str):
    os.makedirs(output_directory, exist_ok=True)
    print(f"Recherche de modèles téléchargeables pour : '{query}'...")

    search_url = 'https://api.sketchfab.com/v3/search'
    params = {
        'type': 'models',
        'q': query,
        'downloadable': True,
        'archives_flavors': True
    }

    try:
        response = requests.get(search_url, params=params, headers=HEADERS)
    except Exception as e:
        print(f"Erreur de connexion lors de la recherche : {e}")
        return

    if response.status_code != 200:
        print(f"Erreur API de recherche : code {response.status_code} - {response.text}")
        return

    results = response.json().get('results', [])
    print(f"{len(results)} spécimen(s) trouvé(s).\n")

    for item in results:
        model_uid = item['uid']
        raw_name = item.get('name', f'model_{model_uid}')
        model_name = sanitize_filename(raw_name)
        folder_name = f"{model_name}_{model_uid}"
        extract_path = os.path.join(output_directory, folder_name)

        print(f"Traitement de : {model_name} (UID: {model_uid})")

        # Obtenir les liens de téléchargement
        download_url = f"https://api.sketchfab.com/v3/models/{model_uid}/download"
        dl_response = requests.get(download_url, headers=HEADERS)

        if dl_response.status_code == 200:
            dl_data = dl_response.json()
            if 'source' in dl_data:
                archive_url = dl_data['source']['url']
                flavor = "source"
            elif 'gltf' in dl_data:
                archive_url = dl_data['gltf']['url']
                flavor = "gltf"
            else:
                print(f"Aucun format adapté trouvé pour {model_uid}.")
                continue

            zip_path = os.path.join(output_directory, f"{folder_name}.zip")

            if os.path.exists(zip_path) and os.path.getsize(zip_path) > 0:
                print(f"L'archive existe déjà localement : {zip_path}")
            else:
                print(f"Téléchargement du modèle ({flavor})...")
                file_response = requests.get(archive_url)
                with open(zip_path, 'wb') as f:
                    f.write(file_response.content)
                print(f"Fichier sauvegardé : {zip_path}")

            # Extraction
            try:
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_path)
                print(f"Archive extraite dans : {extract_path}")
                # Décompression des éventuels zip imbriqués
                extract_nested_zips(extract_path)
            except zipfile.BadZipFile:
                print("L'archive téléchargée est corrompue.")
        else:
            print(f"Erreur d'accès au téléchargement pour {model_uid} (code {dl_response.status_code}).")

    # Recherche et recensement des STL
    find_and_expose_stl_files(output_directory, output_directory)

if __name__ == '__main__':
    search_and_download_specimens(SEARCH_QUERY, OUTPUT_DIR)