import subprocess
import time
import datetime
import os

def catch_date_hour():
    date_hour = datetime.datetime.now()
    return date_hour

def terminal_command(command, input_data=None):
    date_hour = catch_date_hour()
    hour = date_hour.strftime("%H:%M:%S")
    status_file = os.path.join(os.getcwd(), 'status.txt')

    with open(status_file, 'a') as file:
        file.write(f"⏰ [{hour}]: Diretório atual: {os.getcwd()}\n")
        file.write(f"⏰ [{hour}]: Comando '{command}' em execução... ⏳\n")

    print(f'Executando: {command}')

    start_time = time.time()
    terminal_process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if input_data:
        stdout, stderr = terminal_process.communicate(input=input_data)
    else:
        stdout, stderr = terminal_process.communicate()
    
    end_time = time.time()
    total_time = end_time - start_time

    with open(status_file, 'r+') as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            if line.strip() != f"⏰ [{hour}]: Comando '{command}' em execução... ⏳":
                file.write(line)
                print(line)
        
        file.truncate()

        hours, remainder = divmod(total_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        file.write(f"⏰ [{hour}]: Comando '{command}' executado! ✅\n")
        file.write(f"⏰ [{hour}]: Em {int(hours)}h : {int(minutes)}m : {int(seconds)}s ⏳\n\n")
        file.write(f"⏰ [{hour}]: Saída padrão (stdout):\n{stdout}\n")
        
        if stderr:
            file.write(f"⏰ [{hour}]: Saída de erro (stderr):\n{stderr}\n\n")

    return stdout, stderr

# Diretório base onde estão as pastas Variante_XXXX_XXXXX_monomer
base_dir = '/home/insilitox/Documentos/joel/monomeros'
os.chdir(base_dir)

# Itera sobre todas as pastas que correspondem ao padrão Variante_*_monomer
for dir_name in os.listdir(base_dir):
    if dir_name.startswith('Variant') and dir_name.endswith('_monomer'):
        variant_dir = os.path.join(base_dir, dir_name)
        if os.path.isdir(variant_dir):
            print(f"\nProcessando pasta: {variant_dir}")
            
            # Muda para o diretório da variante
            os.chdir(variant_dir)

            # Inicializa o arquivo status.txt na pasta da variante
            status_file = os.path.join(variant_dir, 'status.txt')
            with open(status_file, 'w') as file:
                pass

            with open(status_file, 'a') as file:
                date_hour = catch_date_hour()
                date_hour_format = date_hour.strftime("%d-%m-%Y %H:%M:%S")
                file.write(f"📅🕔 {date_hour_format} 📅🕔\n\n")
                file.write(f"⏳ Por favor, não desligar o PC nem o ar condicionado !!! 😊 \n\n\n")
                file.write(f"📂 Iniciando processamento da pasta: {dir_name}\n")

            # Executa os comandos GROMACS

            terminal_command(f'gmx trjconv -s EM.tpr -f EM.gro -o EM_{dir_name}.pdb', input_data='1\n')

            # terminal_command(f'gmx pdb2gmx -f {dir_name}.pdb -ignh', input_data='1\n1\n')
            # terminal_command('gmx editconf -f conf.gro -d 3.0 -bt dodecahedron -o box.gro')
            # terminal_command('gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o box_solv.gro')
            # terminal_command('gmx grompp -f ions.mdp -c box_solv.gro -maxwarn 2 -p topol.top -o ION.tpr')
            # terminal_command('gmx genion -s ION.tpr -p topol.top -conc 0.15 -neutral -o box_solv_ion.gro', input_data='13\n')
            # terminal_command('gmx grompp -f EM.mdp -c box_solv_ion.gro -maxwarn 2 -p topol.top -o EM.tpr')
            # terminal_command('gmx mdrun -v -deffnm EM')

            # Volta ao diretório base para a próxima iteração
            os.chdir(base_dir)