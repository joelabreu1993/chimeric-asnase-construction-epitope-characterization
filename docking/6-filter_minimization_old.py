import os
import re

# Caminho da pasta onde estão os EM.log
base_dir = "./Variante1337_22442_monomer"

# Regex para capturar o valor de Maximum force
regex_force = re.compile(r"Maximum force\s*=\s*([\d\.Ee\+\-]+)")

# Dicionário para armazenar {arquivo: valor}
forces = {}

# Itera sobre todos os arquivos que começam com EM e terminam com .log
for file in os.listdir(base_dir):
    if file.startswith("EM") and file.endswith(".log"):
        with open(os.path.join(base_dir, file), "r") as f:
            content = f.read()
            match = regex_force.search(content)
            if match:
                forces[file] = float(match.group(1))

# Se não achou nada, encerra
if not forces:
    print("Nenhum arquivo EM*.log encontrado com Maximum force.")
    exit()

# Ordena pelos valores de força
sorted_forces = sorted(forces.items(), key=lambda x: x[1])

# Escreve resumo em arquivo .txt
summary_file = os.path.join(base_dir, "forces_summary.txt")
with open(summary_file, "w") as out:
    out.write("Resumo das forças máximas:\n\n")
    for file, value in sorted_forces:
        out.write(f"{file}: {value:.4f} kJ/mol/nm\n")
    out.write("\nMelhor replicata: "
              f"{sorted_forces[0][0]} com Maximum force = {sorted_forces[0][1]:.4f} kJ/mol/nm\n")

print(f"Resumo salvo em {summary_file}")

# Mantém apenas o melhor arquivo
best_file = sorted_forces[0][0]
for file in forces.keys():
    if file != best_file:
        os.remove(os.path.join(base_dir, file))
        print(f"Arquivo {file} removido.")

print(f"Melhor replicata mantida: {best_file}")
