import requests
import pandas as pd
import time
import os
from Bio import SeqIO

def sanitize_filename(filename):
    return filename.replace("*", "_").replace(":", "_")

def predict_epitopes_lotes(input_file, output_dir_txt, output_dir_excel, batch_size=200):
    api_url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
    alleles = [
        "HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01",
        "HLA-DRB1*07:01", "HLA-DRB1*08:01", "HLA-DRB1*11:01",
        "HLA-DRB1*13:01", "HLA-DRB1*15:01"
    ]
    peptide_length = 15
    method = "netmhciipan_el-4.1"

    # Ler todas as sequências do arquivo FASTA
    try:
        sequences = list(SeqIO.parse(input_file, "fasta"))
    except Exception as e:
        print(f"Erro ao ler o arquivo de entrada: {e}")
        return

    # Criar diretórios de saída se não existirem
    if not os.path.exists(output_dir_txt):
        os.makedirs(output_dir_txt)
    if not os.path.exists(output_dir_excel):
        os.makedirs(output_dir_excel)

    total_sequences = len(sequences)
    processed_sequences = 0
    batch_data = []
    batch_start_id = ""

    for i, sequence_record in enumerate(sequences):
        sequence_id = sanitize_filename(sequence_record.id)
        if i % batch_size == 0:
            batch_start_id = sequence_id  # Nome da primeira sequência do lote

        sequence = str(sequence_record.seq)
        sequence_results = []

        # Para cada alelo, faz a requisição à API (com retry)
        for allele in alleles:
            while True:
                try:
                    payload = {
                        "method": method,
                        "sequence_text": sequence,
                        "allele": allele,
                        "length": str(peptide_length),
                    }
                    response = requests.post(api_url, data=payload)
                    if response.status_code == 200:
                        data = response.text
                        results = []
                        if data.startswith("allele"):
                            for line in data.splitlines():
                                if not line.startswith("allele") and line.strip():
                                    split_line = line.split("\t")
                                    if len(split_line) >= 9:
                                        try:
                                            results.append({
                                                "sequence_id": sequence_id,
                                                "allele": split_line[0],
                                                "seq_num": int(split_line[1]),
                                                "start": int(split_line[2]),
                                                "end": int(split_line[3]),
                                                "length": int(split_line[4]),
                                                "core_peptide": split_line[5],
                                                "peptide": split_line[6],
                                                "score": float(split_line[7]),
                                                "percentile_rank": float(split_line[8]),
                                            })
                                        except ValueError:
                                            continue
                        sequence_results.extend(results)
                        break  # Sai do loop while se bem-sucedido
                    else:
                        print(f"Erro {response.status_code} para o alelo {allele}, tentando novamente...")
                        time.sleep(5)
                except Exception as e:
                    print(f"Erro ao processar o alelo {allele}, tentando novamente: {e}")
                    time.sleep(5)

        # Salvamento dos epítopos em TXT
        if sequence_results:
            try:
                df_sequence = pd.DataFrame(sequence_results)
                output_txt = os.path.join(output_dir_txt, f"{sequence_id}_predicted_epitopes.txt")
                df_sequence.to_csv(output_txt, sep="\t", index=False)
            except Exception as e:
                print(f"Erro ao salvar TXT da sequência {sequence_id}: {e}")

            total_epitopes = len(df_sequence)

            # Calcular densidades e preparar para salvar em Excel
            for allele in alleles:
                subset = df_sequence[df_sequence["allele"] == allele]

                ni_total = subset[subset["percentile_rank"] <= 10].shape[0]
                fi_total = ni_total / total_epitopes if total_epitopes > 0 else 0

                ni_strong = subset[subset["percentile_rank"] <= 2].shape[0]
                fi_strong = ni_strong / total_epitopes if total_epitopes > 0 else 0

                ni_weak = subset[(subset["percentile_rank"] > 2) & (subset["percentile_rank"] <= 10)].shape[0]
                fi_weak = ni_weak / total_epitopes if total_epitopes > 0 else 0

                batch_data.append({
                    "sequence_id": sequence_id,
                    "allele": allele,
                    "ni_total": ni_total,
                    "fi_total": round(fi_total, 4),
                    "ni_strong": ni_strong,
                    "fi_strong": round(fi_strong, 4),
                    "ni_weak": ni_weak,
                    "fi_weak": round(fi_weak, 4),
                    "N": total_epitopes
                })

        processed_sequences += 1
        print(f"Processando sequência {processed_sequences}/{total_sequences}")

        # Ao terminar um lote ou o arquivo todo, salva o Excel do lote
        if (i + 1) % batch_size == 0 or (i + 1) == total_sequences:
            batch_end_id = sequence_id  # Nome da última sequência do lote
            batch_filename = f"{batch_start_id}_to_{batch_end_id}.xlsx"
            batch_filepath = os.path.join(output_dir_excel, batch_filename)
            df_batch = pd.DataFrame(batch_data)
            try:
                df_batch.to_excel(batch_filepath, index=False)
                print(f"Lote salvo em {batch_filepath}")
            except Exception as e:
                print(f"Erro ao salvar o lote Excel: {e}")
            batch_data = []  # Resetar lote

    print("ANÁLISE FINALIZADA COM SUCESSO!")

# Exemplo de uso:
input_file = r'C:/Users/joelabreu1993/OneDrive/Ciências Farmacêuticas - UnB/Projeto de L-asparaginase/Engenharia de proteínas/Trabalho - UFRO/Artigo/Fixed 15/Predição - Percentuais de 10 e 2/Quimeras e originais.fasta'
output_dir_txt = r'C:/Users/joelabreu1993/OneDrive/Ciências Farmacêuticas - UnB/Projeto de L-asparaginase/Engenharia de proteínas/Trabalho - UFRO/Artigo/Fixed 15/Predição - Percentuais de 10 e 2/Quimeras/Resultados_TXT'
output_dir_excel = r'C:/Users/joelabreu1993/OneDrive/Ciências Farmacêuticas - UnB/Projeto de L-asparaginase/Engenharia de proteínas/Trabalho - UFRO/Artigo/Fixed 15/Predição - Percentuais de 10 e 2/Quimeras/Densidade_de_Epitopos_Lotes'

predict_epitopes_lotes(input_file, output_dir_txt, output_dir_excel, batch_size=200)