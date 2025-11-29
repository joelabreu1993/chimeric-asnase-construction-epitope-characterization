import os
import re

# Diretório base onde estão as pastas VarianteXXXX_XXXXX_monomer
base_dir = "."
total = 0
rerun = 0

# Regex para capturar Maximum force
regex_force = re.compile(r"Maximum force\s*=\s*([\d\.Ee\+\-]+)")

# Função para verificar convergência e pegar Maximum force
def verificar_em_log(log_path):
    """Verifica se a minimização convergiu bem pelo conteúdo do EM.log e retorna Fmax."""
    if not os.path.exists(log_path):
        return None, None

    with open(log_path, "r") as f:
        linhas = f.readlines()

    fmax = None
    status = None
    for linha in reversed(linhas):  # lê de baixo para cima
        if "Maximum force" in linha:
            m = regex_force.search(linha)
            if m:
                fmax = float(m.group(1))
        if "Steepest Descents converged to Fmax < 100" in linha:
            status = True
        elif "did not reach the requested Fmax < 100" in linha:
            status = False

    return status, fmax

# Armazena todas as informações para posterior ordenação
info_fmax = []  # lista de tuplas: (folder, fmax, status)

# Loop pelas pastas
for folder in sorted(os.listdir(base_dir)):
    folder_path = os.path.join(base_dir, folder)
    if folder.startswith("Variante") and os.path.isdir(folder_path):
        log_path = os.path.join(folder_path, "EM.log")
        total += 1

        status, fmax = verificar_em_log(log_path)

        if status is True:
            print(f"[OK] {folder}: Convergência atingida (Fmax < 100).")
        elif status is False:
            print(f"[ERRO] {folder}: Não convergiu. Necessita rerun.")
            rerun += 1
        else:
            print(f"[ATENÇÃO] {folder}: Não foi possível interpretar o EM.log.")

        info_fmax.append((folder, fmax, status))

# Escreve o arquivo EM_values.txt
em_values_file = os.path.join(base_dir, "EM_values.txt")
with open(em_values_file, "w") as em_out:
    em_out.write(f"{'Pasta':40s} {'Maximum force (kJ/mol/nm)':>25s}\n")
    em_out.write("-"*70 + "\n")

    # 1) Primeiro as pastas com Fmax > 100 (não convergidas)
    for folder, fmax, status in info_fmax:
        if fmax is not None and fmax > 100:
            em_out.write(f"{folder:40s} {fmax:>25.4f}\n")

    # 2) Depois as pastas com Fmax <= 100 (convergidas)
    for folder, fmax, status in info_fmax:
        if fmax is not None and fmax <= 100:
            em_out.write(f"{folder:40s} {fmax:>25.4f}\n")

    # 3) Pastas sem valor de Fmax
    for folder, fmax, status in info_fmax:
        if fmax is None:
            em_out.write(f"{folder:40s} {'N/A':>25s}\n")

print(f"Ao todo foram verificados {total} pastas.")
print(f"{rerun} pastas precisam de rerun.")
print(f"✅ Valores salvos em {em_values_file}")
