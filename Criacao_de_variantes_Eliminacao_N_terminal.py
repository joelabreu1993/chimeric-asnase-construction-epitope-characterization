from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from itertools import product

# Sequência inicial de aminoácidos
sequence = "SPVHKPQSSSTKYVYTNSNGLNFTQMNPGLPNITIFGTGGTIAGSGSSSTATTGYTAGAVGILDLIDAVPSMLNVSNIAGVQVANVGSEDITSDILISLSKSINKLVCDDSTMAGAVITHGTDTLEETAFFLDATINCGKPVVIVGAMRPSTATSADGPFNLLEAVTVAASPKAVNRGAMVVMNDRIASAYYVTKTNANTMDTFKAIEMGFLGEMISNTPFFFYPPVKPTGKVEFDITATKEIPRVDILYAYEDMHNDTLYSAVENGAQGIVIAGAGAGGVSTSFNHAIEDVINRFKIPVVQSMRTVNGEVPLSDVNSTSAIHIASGYLNPQKSRILLGLLLSEARTLTDIRSVFSLGTVS"

# Definição dos blocos e variantes
mutation_blocks = [
    [("VYTNSNGLNFTQ", 13, 25, "__1."), ("LFSSRRGHLFLQ", 13, 25, "__2."), ("MVKNEDSLVFVQ", 13, 25, "__3."),
     ("YYKEHEGAIY", 13, 23, "__4."), ("LRKLRKRLLR", 13, 23, "__5."), ("AVTLDAGLVY", 13, 23, "__6."), 
     ("QLTTGNGLFL", 13, 23, "__7."), ("LVSLQSGYL", 13, 22, "__8."), ("SFTVG", 13, 18, "__9.")],
    [("TSDILISLSKSINKLVC", 91, 108, "1."), ("KEEFPFALGVQTLPQTC", 91, 108, "2."), ("TEEAPLKLSKAVHKAVL", 91, 108, "3."),
     ("PDATPTELAKLVNKHS", 91, 107, "4."), ("EEYANCHLARAP", 91, 103, "5."), ("KKTPKSPVGVQP", 91, 103, "6."),
     ("KDDPDAPLQPVT", 91, 103, "7."), ("LESFKVSFL", 91, 100, "8.")],
    [("KPTGKVEFDITAT", 227, 240, "1."), ("SPQQVFSTEFEVK", 227, 240, "2."), ("KSGGRTEHPFTVE", 227, 240, "3."),
     ("QCTAEISLADLAT", 227, 240, "4."), ("KDLGEENFKALVL", 227, 240, "5."), ("GYSFTVGGS", 231, 240, "6."), 
     ("NSSTQFEVK", 231, 240, "7.")],
    [("HNDTLYSAVENG",255, 267, "1."), ("GVDYV", 262, 267, "2."), ("VAEYG", 262, 267, "3."),
     ("RVEYG", 262, 267, "4.")],
    [("STSFNHAIEDVINRFK", 281, 297, "1."), ("CEADDGCPKPP", 286, 297, "2."), ("FLGDRDFNQFS", 286, 297, "3."),
     ("FVESKDVCKNY", 286, 297, "4."), ("QDNEDCINRHN", 286, 297, "5.")]
]
    

# Função para gerar mutantes combinando os blocos
def generate_variants(sequence, mutation_blocks):
    all_variants = []
    combinations = product(*mutation_blocks)
    variant_number = 1

    for combination in combinations:
        mutant_seq = sequence
        variant_id = f"L-ASNase_Pc_Variante-{variant_number}."

        for mutation in combination:
            variant, start, end, suffix = mutation
            mutant_seq = mutant_seq[:start] + variant + mutant_seq[end:]
            variant_id += suffix

        mutant_record = SeqRecord(
            Seq(mutant_seq),
            id=variant_id,
            description="_"
        )
        all_variants.append(mutant_record)
        variant_number += 1

    return all_variants

# Função para criar variantes eliminando os 29 primeiros resíduos e ajustando o nome
def generate_variants_with_deletion(original_variants):
    deleted_variants = []

    for record in original_variants:
        mutated_seq = str(record.seq)[29:]  # Eliminar os 29 primeiros resíduos
        # Substituir o primeiro algarismo do sufixo por 0
        variant_id_parts = record.id.split("__")
        variant_id = variant_id_parts[0] + "__0" + variant_id_parts[1][1:]

        deleted_record = SeqRecord(
            Seq(mutated_seq),
            id=variant_id,
            description="_"
        )
        deleted_variants.append(deleted_record)

    return deleted_variants

# Função para eliminar duplicatas
def remove_duplicates(records):
    unique_records = []
    seen_sequences = set()

    for record in records:
        seq_str = str(record.seq)
        if seq_str not in seen_sequences:
            unique_records.append(record)
            seen_sequences.add(seq_str)

    return unique_records

# Gerar as variantes
mutant_records = generate_variants(sequence, mutation_blocks)

# Gerar as variantes com os 29 primeiros resíduos eliminados
deleted_mutant_records = generate_variants_with_deletion(mutant_records)

# Combinar todas as variantes
all_mutant_records = mutant_records + deleted_mutant_records

# Remover duplicatas
unique_mutant_records = remove_duplicates(all_mutant_records)

# Salvar os registros em um arquivo FASTA
output_file = "Quimeras_L-ASNase_PcII.fasta"
SeqIO.write(unique_mutant_records, output_file, "fasta")

print(f"Arquivo FASTA com mutantes combinados gerado: {output_file}")
