import os
import shutil
import subprocess

def run_gold_batch_setup():
    """
    Copia e configura o gold.conf para todas as pastas de variante, executa protonação com gold_utils
    e atualiza as linhas cavity_file e protein_datafile.
    """
    # Encontra o diretório pai (onde estão as pastas 'monomeros' e o 'gold.conf')
    current_directory = os.getcwd()
    monomeros_dir = os.path.join(current_directory, 'monomeros')

    # Define os arquivos a serem copiados
    source_conf = 'gold.conf'
    
    # Caminho do gold_utils
    gold_utils_path = "/home/insilitox/CCDC/ccdc-software/gold/GOLD/gold_utils"

    # Itera sobre os subdiretórios dentro da pasta 'monomeros'
    for variant_name in os.listdir(monomeros_dir):
        # Verifica se o nome da pasta corresponde ao padrão
        if variant_name.startswith('Variante') and variant_name.endswith('_monomer'):
            # Extrai o identificador da variante (ex: "1322_22412")
            variant_id = variant_name.split("Variante")[1].split("_monomer")[0]
            
            # Constrói o caminho completo para a pasta da variante
            variant_folder = os.path.join(monomeros_dir, variant_name)
            
            # Define o nome do arquivo de centroide e do arquivo de destino
            centroid_file_name = 'gold_activesite_aas.txt'
            destination_conf = os.path.join(variant_folder, source_conf)

            # 1. Copia gold.conf para a pasta da variante
            try:
                shutil.copyfile(os.path.join(current_directory, source_conf), destination_conf)
                print(f"gold.conf copiado para: {variant_folder}")
            except FileNotFoundError:
                print(f"Erro: '{source_conf}' não encontrado em: {current_directory}")
                continue

            # 2. Executa a protonação com gold_utils
            input_pdb = f"EM_Variante{variant_id}_monomer.pdb"
            output_pdb = f"EM_Variante{variant_id}_monomer_H.pdb"
            input_pdb_path = os.path.join(variant_folder, input_pdb)
            output_pdb_path = os.path.join(variant_folder, output_pdb)
            
            if not os.path.exists(input_pdb_path):
                print(f"Erro: Arquivo de entrada '{input_pdb}' não encontrado em: {variant_folder}")
                continue
            
            # Adicionando a verificação e remoção do arquivo de saída se ele já existir
            if os.path.exists(output_pdb_path):
                print(f"Aviso: O arquivo de saída '{output_pdb_path}' já existe. Removendo...")
                os.remove(output_pdb_path)

            # Configura o comando gold_utils
            gold_utils_cmd = [
                gold_utils_path,
                "-protonate",
                "-i", input_pdb_path,
                "-o", output_pdb_path
            ]
            
            try:
                # Executa o comando usando o ambiente padrão
                result = subprocess.run(gold_utils_cmd, capture_output=True, text=True, check=True)
                print(f"Protonação concluída para {variant_name}: {output_pdb_path}")
                print(f"Saída do gold_utils: {result.stdout}")
                if result.stderr:
                    print(f"Avisos/erros do gold_utils: {result.stderr}")
            except subprocess.CalledProcessError as e:
                print(f"Erro ao executar gold_utils para {variant_name}: {e}")
                print(f"Saída de erro: {e.stderr}")
                continue
            except FileNotFoundError:
                print(f"Erro: '{gold_utils_path}' não encontrado. Verifique o caminho do gold_utils.")
                continue

            # 3. Verifica a existência de gold_activesite_aas.txt
            centroid_file_path = os.path.join(variant_folder, centroid_file_name)
            if not os.path.exists(centroid_file_path):
                print(f"Erro: '{centroid_file_name}' não encontrado em: {variant_folder}")
                continue

            # 4. Altera o gold.conf com as novas informações
            try:
                with open(destination_conf, 'r') as f:
                    conf_lines = f.readlines()

                new_conf_lines = []
                protein_datafile_updated = False
                cavity_file_updated = False
                for line in conf_lines:
                    if line.strip().startswith('cavity_file ='):
                        new_conf_lines.append(f"cavity_file = {centroid_file_name}\n")
                        cavity_file_updated = True
                    elif line.strip().startswith('protein_datafile ='):
                        new_conf_lines.append(f"protein_datafile = {output_pdb_path}\n")
                        protein_datafile_updated = True
                    else:
                        new_conf_lines.append(line)

                if not protein_datafile_updated:
                    print(f"Aviso: Linha 'protein_datafile =' não encontrada em gold.conf para {variant_name}")
                if not cavity_file_updated:
                    print(f"Aviso: Linha 'cavity_file =' não encontrada em gold.conf para {variant_name}")

                with open(destination_conf, 'w') as f:
                    f.writelines(new_conf_lines)
                
                print(f"gold.conf em '{variant_folder}' atualizado com sucesso!")
                print("-" * 20)
            
            except Exception as e:
                print(f"Erro ao processar a pasta '{variant_name}': {e}")
                
    print("Processo de configuração finalizado para todas as variantes.")

# Chama a função para executar a tarefa
run_gold_batch_setup()