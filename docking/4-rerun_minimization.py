import os
import subprocess

# Diretório base onde estão as pastas VarianteXXXX_XXXXX_monomer
base_dir = "."
total = 0
rerun = 0

def verificar_em_log(log_path):
    """Verifica se a minimização convergiu bem pelo conteúdo do EM.log."""
    with open(log_path, "r") as f:
        linhas = f.readlines()

    for linha in reversed(linhas):  # lê de baixo para cima
        if "Steepest Descents converged to Fmax < 100" in linha:
            return True
        if "did not reach the requested Fmax < 100" in linha:
            return False
    return None  # caso não encontre nada

# Loop pelas pastas
for folder in os.listdir(base_dir):
    if folder.startswith("Variante") and os.path.isdir(os.path.join(base_dir, folder)):
        log_path = os.path.join(base_dir, folder, "EM.log")
        total += 1

        if not os.path.exists(log_path):
            print(f"[AVISO] {folder} não contém EM.log")
            continue

        status = verificar_em_log(log_path)

        if status is True:
            print(f"[OK] {folder}: Convergência atingida (Fmax < 100).")
        elif status is False:
            # print(f"[ERRO] {folder}: Não convergiu. Reexecutando minimização...")
            # subprocess.run("gmx grompp -f EM.mdp -c box_solv_ion.gro -maxwarn 2 -p topol.top -o EM_1.tpr", cwd=os.path.join(base_dir, folder), shell=True)
            # subprocess.run("gmx mdrun -v -deffnm EM_1", cwd=os.path.join(base_dir, folder), shell=True)

            # subprocess.run("gmx grompp -f EM.mdp -c box_solv_ion.gro -maxwarn 2 -p topol.top -o EM_2.tpr", cwd=os.path.join(base_dir, folder), shell=True)
            # subprocess.run("gmx mdrun -v -deffnm EM_2", cwd=os.path.join(base_dir, folder), shell=True)

            # subprocess.run("gmx grompp -f EM.mdp -c box_solv_ion.gro -maxwarn 2 -p topol.top -o EM_3.tpr", cwd=os.path.join(base_dir, folder), shell=True)
            # subprocess.run("gmx mdrun -v -deffnm EM_3", cwd=os.path.join(base_dir, folder), shell=True)

            # subprocess.run("gmx grompp -f EM.mdp -c box_solv_ion.gro -maxwarn 2 -p topol.top -o EM_4.tpr", cwd=os.path.join(base_dir, folder), shell=True)
            # subprocess.run("gmx mdrun -v -deffnm EM_4", cwd=os.path.join(base_dir, folder), shell=True)

            # subprocess.run("gmx grompp -f EM.mdp -c box_solv_ion.gro -maxwarn 2 -p topol.top -o EM_5.tpr", cwd=os.path.join(base_dir, folder), shell=True)
            # subprocess.run("gmx mdrun -v -deffnm EM_5", cwd=os.path.join(base_dir, folder), shell=True)
            rerun += 1
        else:
            print(f"[ATENÇÃO] {folder}: Não foi possível interpretar o EM.log.")

print(f'Ao todo foram verificados {total} arquivos.')
print(f'{rerun} arquivos foram recalculados.')

# Antes de recalcular:
# Ao todo foram verificados 112 arquivos.
# 71 arquivos foram recalculados.

# Depois de recalcular:
# Ao todo foram verificados 112 arquivos.
# 21 arquivos foram recalculados.

