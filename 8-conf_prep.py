import os
import shutil

def setup_gold_for_variant():
    """
    Copia o gold.conf, atualiza as coordenadas e o nome da proteína
    para uma única pasta de variante, definida diretamente no código.
    """
    # Define o nome da pasta da variante a ser processada
    # **CORRIGIDO**: Removido o "./monomeros/" do caminho.
    variant_name = "Variante1322_22412_monomer"
    
    # Define os caminhos dos arquivos
    current_directory = os.getcwd()
    source_conf = 'gold.conf'
    centroid_file_name = 'normal_centroid.txt'
    
    # **CORRIGIDO**: O caminho da pasta variante agora é construído corretamente.
    variant_folder = os.path.join(current_directory, 'monomeros', variant_name)
    
    # Verifica se a pasta existe
    if not os.path.exists(variant_folder):
        print(f"Erro: A pasta '{variant_folder}' não foi encontrada.")
        return

    # 1. Copia gold.conf para a pasta da variante
    destination_conf = os.path.join(variant_folder, source_conf)
    # **CORRIGIDO**: Caminho de origem do gold.conf agora é o diretório atual do script.
    shutil.copyfile(os.path.join(current_directory, source_conf), destination_conf)
    print(f"gold.conf copiado para: {variant_folder}")

    # 2. Lê as coordenadas do arquivo normal_centroid.txt
    centroid_file_path = os.path.join(variant_folder, centroid_file_name)
    if os.path.exists(centroid_file_path):
        with open(centroid_file_path, 'r') as f:
            lines = f.readlines()
            # A terceira linha tem as coordenadas (índice 2)
            coordinates_line = lines[1].strip()
            # **CORRIGIDO**: Use split() sem argumentos para lidar com múltiplos espaços
            coordinates = ' '.join(coordinates_line.split())
    else:
        print(f"Erro: '{centroid_file_name}' não encontrado em: {variant_folder}")
        return

    # 3. Altera o gold.conf com as novas informações
    with open(destination_conf, 'r') as f:
        conf_content = f.read()

    # Encontra o nome do arquivo da proteína
    protein_file_name = f"{variant_name}.pdb"
    
    # **CORRIGIDO**: A linha de substituição agora usa o caminho completo correto.
    # O caminho completo é a pasta da variante + o nome do arquivo da proteína.
    
    # Encontra a linha de 'origin' e a substitui com as novas coordenadas
    new_conf_lines = []
    for line in conf_content.splitlines():
        if line.strip().startswith('origin ='):
            new_conf_lines.append(f"origin = {coordinates}")
        elif line.strip().startswith('protein_datafile ='):
            # Encontra a linha 'protein_datafile' e a substitui com o novo caminho
            new_conf_lines.append(f"protein_datafile = {os.path.join(variant_folder, protein_file_name)}")
        else:
            new_conf_lines.append(line)
            
    conf_content = "\n".join(new_conf_lines)

    # Salva o arquivo gold.conf atualizado
    with open(destination_conf, 'w') as f:
        f.write(conf_content)
    
    print(f"gold.conf em '{variant_folder}' atualizado com sucesso!")

# Chama a função para executar a tarefa
setup_gold_for_variant()
