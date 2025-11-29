import os
import shutil

# Define o diretório atual como a pasta 'monomeros'
monomeros_dir = os.getcwd()

# Lista de arquivos e pastas a serem copiados
items_to_copy = ['charmm36-jul2022.ff', 'EM.mdp', 'ions.mdp']

# Itera sobre todas as pastas que correspondem ao padrão Variante_*_monomer
for dir_name in os.listdir(monomeros_dir):
    if dir_name.startswith('Variante') and dir_name.endswith('_monomer'):
        dest_dir = os.path.join(monomeros_dir, dir_name)
        # Verifica se é uma pasta
        if os.path.isdir(dest_dir):
            for item in items_to_copy:
                source_path = os.path.join(monomeros_dir, item)
                dest_path = os.path.join(dest_dir, item)
                # Copia pastas (recursivamente) ou arquivos
                if os.path.isdir(source_path):
                    shutil.copytree(source_path, dest_path, dirs_exist_ok=True)
                else:
                    shutil.copy2(source_path, dest_path)
                print(f'Copiado {item} para {dest_dir}')