"""
Módulo para a Função 1c: Enviar (preparar) sequências para modelagem.
"""

import os
import re
from Bio import SeqIO

def enviar_para_modelagem(dir_leitura_fasta, dir_escrita_model):
    """
    Lê TODOS os arquivos FASTA de 'Funcao1_Filtrar', 'Funcao2a_Separar' ou 'input' 
    e salva os FASTAs individuais em subpastas dentro de 'Funcao2b_FastasIndividuais'.
    """
    try:
        fasta_files = [f for f in os.listdir(dir_leitura_fasta) if f.endswith(".fasta")]
        if not fasta_files:
            print(f"\n[Atenção] Nenhum arquivo FASTA encontrado em '{dir_leitura_fasta}'. Etapa 1c pulada.")
            return None

        print(f"Encontrados {len(fasta_files)} arquivos .fasta para processar.")

        pasta_modelagem = dir_escrita_model
        resultados_modelagem = {}

        for arquivo_fasta in fasta_files:
            print(f"\n  --- Processando arquivo: {arquivo_fasta} ---")
            caminho_fasta = os.path.join(dir_leitura_fasta, arquivo_fasta)

            seqs = list(SeqIO.parse(caminho_fasta, "fasta"))
            if not seqs:
                print(f"  Atenção: o arquivo {arquivo_fasta} não contém sequências. Pulando.")
                continue

            nome_base = os.path.splitext(arquivo_fasta)[0]
            subpasta = os.path.join(pasta_modelagem, nome_base)
            os.makedirs(subpasta, exist_ok=True)

            arquivos_individuais = []

            for i, record in enumerate(seqs, start=1):
                nome_seq = re.sub(r'[^A-Za-z0-9_-]', '_', record.id)
                arquivo_individual = os.path.join(subpasta, f"{nome_seq}.fasta")

                with open(arquivo_individual, 'w') as f:
                    SeqIO.write(record, f, "fasta")

                arquivos_individuais.append(arquivo_individual)

            print(f"  {len(arquivos_individuais)} arquivos individuais gerados em: {subpasta}")
            resultados_modelagem[arquivo_fasta] = arquivos_individuais

        print("\nPreparação para modelagem (1c) concluída!")
        return resultados_modelagem

    except Exception as e:
        print(f"\nErro ao preparar sequências para modelagem: {e}")
        return None