import os
import subprocess
import time
import datetime

def update_status_header(file_path, status):
    """
    Atualiza a mensagem inicial do arquivo de status.
    'running' para iniciar, 'completed' para finalizar.
    """
    with open(file_path, 'r+') as file:
        lines = file.readlines()
        file.seek(0)
        
        # Mantém a primeira linha (data/hora)
        file.write(lines[0])
        file.write("\n")
        
        # Escreve a nova mensagem de status
        if status == 'running':
            file.write("⏳ Por favor, não desligar o PC nem o ar condicionado !!! 😊\n\n")
        elif status == 'completed':
            file.write("✅ Todas as tarefas foram finalizadas. O PC pode ser desligado. 😊\n\n")
            
        # Adiciona o restante das linhas do log
        file.writelines(lines[3:])
        file.truncate()

def run_all_variants():
    """
    Navega por cada pasta de variante e executa o docking em cada uma.
    """
    # Diretório base onde estão as pastas Variante_XXXX_XXXXX_monomer
    base_dir = '/home/insilitox/Documentos/joel/monomeros'
    
    # Define o nome do arquivo de status global
    global_status_file = os.path.join(base_dir, 'docking_status.txt')
    
    # Use o caminho completo para o executável 'gold_auto'
    docking_command = ['/home/insilitox/CCDC/ccdc-software/gold/GOLD/bin/gold_auto', 'gold.conf']

    # Mensagem de início para o status global
    start_time_str = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    
    with open(global_status_file, 'w') as file:
        file.write(f"📅🕔 {start_time_str} 📅🕔\n\n")
        file.write("⏳ Por favor, não desligar o PC nem o ar condicionado !!! 😊\n\n")

    # Itera sobre todas as pastas que correspondem ao padrão
    for dir_name in os.listdir(base_dir):
        if dir_name.startswith('Variante') and dir_name.endswith('_monomer'):
            variant_dir = os.path.join(base_dir, dir_name)
            
            if os.path.isdir(variant_dir):
                print(f"\nProcessando pasta: {variant_dir}")
                
                # Registra o início do processamento no arquivo de status global
                with open(global_status_file, 'a') as file:
                    file.write(f"📂 Iniciando processamento da pasta: {dir_name}\n")
                
                try:
                    # Executa o comando de docking do GOLD
                    result = subprocess.run(
                        docking_command,
                        cwd=variant_dir,
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    
                    # Se o comando rodou sem erros
                    with open(global_status_file, 'a') as file:
                        file.write("    ✅ Docking concluído com sucesso.\n\n")

                except FileNotFoundError:
                    with open(global_status_file, 'a') as file:
                        file.write("❌ ERRO: Comando 'gold_auto' não encontrado.\n")
                        file.write("Verifique se ele está no seu $PATH ou use o caminho completo.\n")
                        file.write("-" * 20 + "\n\n")
                except subprocess.CalledProcessError as e:
                    with open(global_status_file, 'a') as file:
                        file.write("❌ ERRO: O comando de docking retornou um erro.\n")
                        file.write(f"Saída de erro:\n{e.stderr}\n")
                        file.write("-" * 20 + "\n\n")
                except Exception as e:
                    with open(global_status_file, 'a') as file:
                        file.write(f"❌ ERRO INESPERADO: {e}\n")
                        file.write("-" * 20 + "\n\n")
            
    # Mensagem final, após todas as tarefas
    update_status_header(global_status_file, 'completed')
        
    print("\nProcesso de docking finalizado para todas as variantes. Verifique o arquivo 'docking_status.txt'.")
    
run_all_variants()
