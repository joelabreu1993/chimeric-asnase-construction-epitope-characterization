import os
from pymol import cmd

# Diretório principal
monomeros_dir = "./monomeros"

# --- Reinicializar o PyMOL para um estado limpo ---
cmd.reinitialize()

# --- Carregar o arquivo PDB comum (na raiz) ---
cmd.load("6v2a_monomer.pdb")

# Listar todas as subpastas em ./monomeros que seguem o padrão "VarianteXXXX_XXXXX_monomer"
pastas = [d for d in os.listdir(monomeros_dir) if os.path.isdir(os.path.join(monomeros_dir, d)) and d.startswith("Variante") and d.endswith("_monomer")]

for pasta in pastas:
    # Extrair o variant do nome da pasta (ex: "1322_22412")
    variant = pasta.split("Variante")[1].split("_monomer")[0]
    
    # Nome do modelo e caminho completo da pasta
    full_pasta = os.path.join(monomeros_dir, pasta)
    model_name = f"EM_Variante{variant}_monomer"
    
    # --- Carregar o arquivo PDB da variante ---
    pdb_path = os.path.join(full_pasta, f"EM_Variante{variant}_monomer.pdb")
    if not os.path.exists(pdb_path):
        print(f"Arquivo PDB não encontrado em {full_pasta}: {pdb_path}")
        continue
    cmd.load(pdb_path, model_name)
    
    # --- Alinhamento ---
    try:
        rmsd = cmd.align("6v2a_monomer", model_name)[0]
        print(f"RMSD do alinhamento para {variant}: {rmsd:.3f}")
    except:
        print(f"Erro no alinhamento para {variant}. Pulando...")
        cmd.delete(model_name)
        continue
    
    # --- Remover moléculas de água ---
    cmd.remove("resn HOH")
    cmd.remove("resn WAT")
    
    # --- Seleção do ligante (ASN 401 chain A) ---
    cmd.select("lig", "hetatm and chain A and resn ASN and resi 401")
    lig_atoms = cmd.count_atoms("lig")
    print(f"Número de átomos na seleção 'lig' para {variant}: {lig_atoms}")
    if lig_atoms == 0:
        print(f"Erro: Ligante ASN 401 (cadeia A) não encontrado para {variant}.")
        cmd.delete(model_name)
        continue
    cmd.show("sticks", "lig")
    
    # --- Seleções para pocket ---
    cmd.select("pocket_6v2a", "byres (polymer within 5.0 of lig)")
    print(f"Seleção pocket_6v2a para {variant}: {cmd.count_atoms('pocket_6v2a')} átomos")
    cmd.select("pocket_model", f"{model_name} within 5.0 of pocket_6v2a")
    print(f"Seleção pocket_model para {variant}: {cmd.count_atoms('pocket_model')} átomos")
    cmd.select("near_model", "byres (polymer and pocket_model)")
    print(f"Seleção near_model para {variant}: {cmd.count_atoms('near_model')} átomos")
    
    # --- Obter átomos CA da seleção ---
    atoms = cmd.get_model("near_model and name CA").atom
    
    if not atoms:
        print(f"Nenhum átomo CA encontrado na seleção 'near_model' para {variant}.")
        cmd.delete(model_name)
        cmd.delete("lig")
        cmd.delete("pocket_6v2a")
        cmd.delete("pocket_model")
        cmd.delete("near_model")
        continue
    
    # --- Calcular centróide ---
    x = sum(a.coord[0] for a in atoms) / len(atoms)
    y = sum(a.coord[1] for a in atoms) / len(atoms)
    z = sum(a.coord[2] for a in atoms) / len(atoms)
    print(f"Centróide (X, Y, Z) para {variant}: ({x}, {y}, {z})")
    
    # --- Criar pseudoátomo no centroide ---
    cmd.pseudoatom("centro_pocket", pos=[x, y, z])
    cmd.show("spheres", "centro_pocket")
    cmd.color("red", "centro_pocket")
    cmd.set("sphere_scale", 1.0, "centro_pocket")
    
    # --- Selecionar resíduos em 8 Å do centroide ---
    cmd.select("pocket_residues", f"(byres ({model_name} within 8 of centro_pocket) and polymer.protein)")
    
    # --- Mostrar proteínas como cartoon e colorir, excluindo o ligante ---
    cmd.show("cartoon", model_name)
    cmd.color("slate", f"{model_name} and not lig")
    cmd.show("cartoon", "6v2a_monomer")
    cmd.color("wheat", "6v2a_monomer and not lig")
    
    # --- Ajustar zoom ---
    cmd.zoom("centro_pocket", buffer=15)
    
    # --- Salvar informações no arquivo coordinates.txt ---
    residues = []
    for a in cmd.get_model("pocket_residues and name CA").atom:
        residues.append((a.resi, a.resn, a.coord[0], a.coord[1], a.coord[2]))
    n_residues = len(residues)
    
    coordinates_path = os.path.join(full_pasta, "coordinates.txt")
    with open(coordinates_path, "w") as f:
        print(f"Escrevendo coordinates.txt em {coordinates_path}")
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
    gold_activesite_path = os.path.join(full_pasta, "gold_activesite_aas.txt")
    with open(gold_activesite_path, "w") as f:
        print(f"Escrevendo gold_activesite_aas.txt em {gold_activesite_path}")
        f.write("> <Gold.Protein.ActiveResidues>\n")
        f.write(" ".join(residues_unique) + "\n")
    
    # --- Ajustes visuais para publicação ---
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_trace_gain", 0.1)
    cmd.set("ray_trace_disco_factor", 1)
    cmd.set("ray_trace_fog", 0.5)
    
    # --- Colorir o ligante por último ---
    cmd.color("yellow", "lig")
    cmd.rebuild()
    cmd.refresh()
    
    # --- Salvar imagem sem ray tracing ---
    image_path = os.path.join(full_pasta, "centro_pocket.png")
    cmd.png(image_path, width=2000, height=1800, dpi=300, ray=0)
    print(f"Imagem 'centro_pocket.png' salva com sucesso para {variant} em {full_pasta}!\n")
    print("="*20)
    
    # --- Limpar o modelo da variante atual para o próximo loop ---
    cmd.delete(model_name)
    cmd.delete("lig")
    cmd.delete("pocket_6v2a")
    cmd.delete("pocket_model")
    cmd.delete("near_model")
    cmd.delete("centro_pocket")
    cmd.delete("pocket_residues")