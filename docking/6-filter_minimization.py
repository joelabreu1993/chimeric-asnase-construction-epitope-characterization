import os
import re

# você está na pasta "monomeros"
base_dir = "."

# regex para capturar o valor de Maximum force nos .log
regex_force = re.compile(r"Maximum force\s*=\s*([\d\.Ee\+\-]+)")

# cria/zera o arquivo global de status
global_status_path = os.path.join(base_dir, "replicated_status.txt")
with open(global_status_path, "w") as g:
    g.write("Resumo geral das replicatas:\n\n")

# percorre todas as pastas Variante*_monomer
for dir_name in sorted(os.listdir(base_dir)):
    variant_dir = os.path.join(base_dir, dir_name)
    if not (os.path.isdir(variant_dir) and dir_name.startswith("Variante") and dir_name.endswith("_monomer")):
        continue

    print(f"\n📂 Processando: {dir_name}")

    # 1) Ler todos os EM*.log e coletar Fmax
    forces = {}  # {base_sem_ext: valor}
    log_files = [f for f in os.listdir(variant_dir) if f.startswith("EM") and f.endswith(".log") and f != "EM.mdp"]

    for log_file in log_files:
        log_path = os.path.join(variant_dir, log_file)
        with open(log_path, "r", errors="ignore") as f:
            content = f.read()
        m = regex_force.search(content)
        if m:
            base_no_ext = os.path.splitext(log_file)[0]  # EM, EM_1, EM_2, ...
            forces[base_no_ext] = float(m.group(1))

    if not forces:
        print("⚠️  Nenhum EM*.log com 'Maximum force' encontrado. Pulando pasta.")
        continue

    # 2) Ordenar por Fmax e salvar resumo local
    sorted_forces = sorted(forces.items(), key=lambda x: x[1])  # [(base, valor), ...]
    best_base, best_f = sorted_forces[0]

    summary_path = os.path.join(variant_dir, "forces_summary.txt")
    with open(summary_path, "w") as out:
        out.write("Resumo das forças máximas (kJ/mol/nm):\n\n")
        for base, val in sorted_forces:
            out.write(f"{base}.log: {val:.4f}\n")
        out.write(f"\nMelhor replicata: {best_base}.log  (Maximum force = {best_f:.4f} kJ/mol/nm)\n")

    print(f"📝 Resumo salvo em: {summary_path}")
    print(f"✅ Melhor replicata: {best_base}.log (Fmax={best_f:.2f})")

    # também salvar no arquivo global
    with open(global_status_path, "a") as g:
        g.write(f"📂 {dir_name}\n")
        g.write(f"   📝 Resumo salvo em: {summary_path}\n")
        g.write(f"   ✅ Melhor replicata: {best_base}.log (Fmax={best_f:.2f})\n\n")

    # 3) Apagar todos os arquivos EM*.* que NÃO pertençam à replicata escolhida
    em_pattern_files = [f for f in os.listdir(variant_dir) if re.match(r"^EM(?:_\d+)?\..+$", f)]
    for fname in em_pattern_files:
        # sempre manter o EM.mdp
        if fname == "EM.mdp":
            continue
        base_no_ext = os.path.splitext(fname)[0]
        if base_no_ext != best_base:
            try:
                os.remove(os.path.join(variant_dir, fname))
                print(f"🗑️  Removido: {fname}")
            except FileNotFoundError:
                pass

    # 4) Renomear os arquivos da replicata escolhida EM_X.* -> EM.*
    if best_base != "EM":
        chosen_files = [f for f in os.listdir(variant_dir) if os.path.splitext(f)[0] == best_base]
        for fname in chosen_files:
            if fname == "EM.mdp":  # não renomear o .mdp
                continue
            old_path = os.path.join(variant_dir, fname)
            ext = os.path.splitext(fname)[1]  # .log, .edr, .gro, .tpr, .trr, ...
            new_path = os.path.join(variant_dir, "EM" + ext)
            if os.path.exists(new_path):
                os.remove(new_path)
            os.rename(old_path, new_path)
            print(f"✏️  Renomeado: {fname}  →  EM{ext}")

print("\n🎉 Pronto! Todas as pastas Variante foram processadas.")
print(f"📑 Resumo geral salvo em: {global_status_path}")
