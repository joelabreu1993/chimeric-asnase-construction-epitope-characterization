from pymol import cmd

# --- Reinicializar o PyMOL para um estado limpo ---
cmd.reinitialize()

# --- Carregar arquivos PDB ---
cmd.load("6v2a_monomer.pdb")
cmd.load("EM_Variante1322_22412_monomer.pdb")

# --- Alinhamento ---
cmd.align("6v2a_monomer", "EM_Variante1322_22412_monomer")

# --- Remover moléculas de água ---
cmd.remove("resn HOH")
cmd.remove("resn WAT")

# --- Seleção do ligante (ASN 401 chain A) ---
cmd.select("lig", "hetatm and chain A and resn ASN and resi 401")
print(f"Número de átomos na seleção 'lig': {cmd.count_atoms('lig')}")  # Depuração
cmd.show("sticks", "lig")  # Mostrar como sticks

# --- Seleções para pocket ---
cmd.select("pocket_6v2a", "byres (polymer within 3.5 of lig)")
cmd.select("pocket_model", "EM_Variante1322_22412_monomer within 3.5 of pocket_6v2a")
cmd.select("near_model", "byres (polymer and pocket_model)")

# --- Obter átomos CA da seleção ---
atoms = cmd.get_model("near_model and name CA").atom

if not atoms:
    print("Nenhum átomo CA encontrado na seleção 'near_model'.")
else:
    # --- Calcular centróide ---
    x = sum(a.coord[0] for a in atoms) / len(atoms)
    y = sum(a.coord[1] for a in atoms) / len(atoms)
    z = sum(a.coord[2] for a in atoms) / len(atoms)
    print(f"Centróide (X, Y, Z): ({x}, {y}, {z})")

    # --- Criar pseudoátomo no centroide ---
    cmd.pseudoatom("centro_pocket", pos=[x, y, z])
    cmd.show("spheres", "centro_pocket")
    cmd.color("red", "centro_pocket")
    cmd.set("sphere_scale", 1.0, "centro_pocket")

    # --- Selecionar resíduos em 10 Å do centroide ---
    cmd.select("pocket_residues", "(byres (EM_Variante1322_22412_monomer within 8 of centro_pocket) and polymer.protein)")

    # --- Mostrar proteínas como cartoon e colorir, excluindo o ligante ---
    cmd.show("cartoon", "EM_Variante1322_22412_monomer")
    cmd.color("slate", "EM_Variante1322_22412_monomer and not lig")
    cmd.show("cartoon", "6v2a_monomer")
    cmd.color("wheat", "6v2a_monomer and not lig")

    # --- Ajustar zoom ---
    cmd.zoom("centro_pocket", buffer=15)

    # --- Salvar informações no arquivo coordinates.txt ---
    residues = []
    for a in cmd.get_model("pocket_residues and name CA").atom:
        residues.append((a.resi, a.resn, a.coord[0], a.coord[1], a.coord[2]))
    n_residues = len(residues)

    with open("coordinates.txt", "w") as f:
        f.write(f"Coordinates X, Y e Z of centroid:\n")
        f.write(f"{x:10.3f} {y:10.3f} {z:10.3f}\n\n")
        f.write(f"{n_residues} aminoácids within a 8 Å radius of the centroid:\n\n")
        f.write(f"{'Resi':>5} {'Resn':>5} {'X':>10} {'Y':>10} {'Z':>10}\n")
        f.write("-" * 45 + "\n")
        for resi, resn, cx, cy, cz in residues:
            f.write(f"{resi:>5} {resn:>5} {cx:10.3f} {cy:10.3f} {cz:10.3f}\n")
        f.write("\n" + "="*50 + "\n\n")

    # --- Preparar lista única de resíduos para Gold ---
    residues_unique = []
    for resi, resn, _, _, _ in residues:
        tag = f"{resn}{resi}"
        if tag not in residues_unique:
            residues_unique.append(tag)

    # --- Salvar informações no arquivo gold_activesite_aas.txt ---
    with open("gold_activesite_aas.txt", "w") as f:
        f.write("> <Gold.Protein.ActiveResidues>\n")
        f.write(" ".join(residues_unique) + "\n")

    # --- Ajustes visuais para publicação ---
    cmd.bg_color("white")          # fundo branco
    cmd.set("ray_trace_mode", 1)   # iluminação realista
    cmd.set("antialias", 2)        # suavizar bordas
    cmd.set("ray_shadows", 0)      # remover sombras
    cmd.set("ray_trace_gain", 0.1)
    cmd.set("ray_trace_disco_factor", 1)
    cmd.set("ray_trace_fog", 0.5)

    # --- Colorir o ligante por último para evitar conflitos ---
    cmd.color("yellow", "lig")
    cmd.rebuild()  # Reconstruir a cena após coloração
    cmd.refresh()  # Forçar atualização visual

    # --- Salvar imagem intermediária sem ray tracing para depuração ---
    cmd.png("centro_pocket.png", width=2000, height=1800, dpi=300, ray=0)
    print("Imagem de depuração 'centro_pocket.png' salva com sucesso!")