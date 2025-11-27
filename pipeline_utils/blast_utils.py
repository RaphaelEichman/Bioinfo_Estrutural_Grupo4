"""
Módulo para a Função 2: Rodar BLASTp.
"""

import os
import subprocess
from Bio import SeqIO
import config

def rodar_blast(dir_leitura_fasta, dir_escrita_blast, automatico=False):
    """
    Lê arquivos FASTA e roda BLASTp.
    Se automatico=True, processa todos os arquivos da pasta sem perguntar.
    """
    try:
        fasta_files = [f for f in os.listdir(dir_leitura_fasta) if f.endswith('.fasta')]
        
        if not fasta_files:
            print(f"\n[Atenção] Nenhum arquivo .fasta encontrado em '{dir_leitura_fasta}'.")
            return

        # LÓGICA DE SELEÇÃO (Pergunta 2)
        if automatico:
            # Pula a pergunta e seleciona tudo
            print(f"-> Processando automaticamente {len(fasta_files)} arquivo(s) em '{os.path.basename(dir_leitura_fasta)}'...")
            arquivos_escolhidos = fasta_files
        else:
            # Modo manual (Faz a pergunta 2)
            print(f"\nArquivos .fasta encontrados em '{os.path.basename(dir_leitura_fasta)}':")
            for f in fasta_files:
                print(f"  - {f}")
                
            print("\nQuais arquivos desta pasta você deseja processar?")
            for i, f in enumerate(fasta_files, start=1):
                print(f"[{i}] {f}")

            selecao = input("\nDigite os números (ex: 1,3) ou 0 para todos: ").strip()

            if selecao == "0":
                arquivos_escolhidos = fasta_files
            else:
                try:
                    indices = [int(x.strip()) for x in selecao.split(',') if x.strip().isdigit()]
                    arquivos_escolhidos = [fasta_files[i - 1] for i in indices if 1 <= i <= len(fasta_files)]
                except Exception:
                    print("Seleção inválida. Pulando BLAST.")
                    return

        if not arquivos_escolhidos:
            print("Nenhum arquivo selecionado. Pulando BLAST.")
            return

        print("\n--- Iniciando BLASTp ---")
        for fasta in arquivos_escolhidos:
            entrada = os.path.join(dir_leitura_fasta, fasta)
            saida = os.path.join(dir_escrita_blast, f"{fasta}_blast.tsv")

            # Verifica se já existe para não refazer (opcional, mas útil)
            if os.path.exists(saida) and os.path.getsize(saida) > 0:
                print(f"  -> {fasta} já processado. Pulando.")
                continue

            print(f"Rodando BLASTp para {fasta}...")

            comando = [
                "blastp",
                "-query", entrada,
                "-db", "pdb",
                "-remote",
                "-evalue", "1e-5",
                "-max_target_seqs", str(config.blast_max_target_seqs), 
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
                "-out", saida
            ]

            subprocess.run(comando)
            print(f"Resultado salvo em: {saida}")

        print(f"\nBLAST concluído para pasta: {os.path.basename(dir_leitura_fasta)}")
        
    except Exception as e:
        print(f"\nErro inesperado no BLAST: {e}")