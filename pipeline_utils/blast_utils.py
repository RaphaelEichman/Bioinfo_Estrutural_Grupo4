# pipeline_utils/blast_utils.py
"""
Módulo para a Função 2: Rodar BLASTp.
(Versão otimizada sem "splitting" de arquivos)
"""

import os
import subprocess
from Bio import SeqIO
import config  # <-- ADICIONADO IMPORT

def rodar_blast(dir_leitura_fasta, dir_escrita_blast):
    """
    Lê arquivos FASTA de 'dir_leitura_fasta' e salva os
    resultados do BLAST em 'dir_escrita_blast'.
    
    Esta versão roda o BLAST para cada ARQUIVO .fasta inteiro.
    """
    try:
        # 1. Encontra os arquivos .fasta na pasta de leitura
        fasta_files = [f for f in os.listdir(dir_leitura_fasta) if f.endswith('.fasta')]
        
        if not fasta_files:
            print(f"\n[Atenção] Nenhum arquivo .fasta encontrado em '{dir_leitura_fasta}'.")
            return

        print(f"\nArquivos .fasta encontrados em '{dir_leitura_fasta}':")
        for f in fasta_files:
            print(f"  - {f}")

        # 2. Pergunta ao usuário quais processar
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

        # 3. Processa os arquivos selecionados
        print("\n--- Iniciando BLASTp ---")
        for fasta in arquivos_escolhidos:
            entrada = os.path.join(dir_leitura_fasta, fasta)
            saida = os.path.join(dir_escrita_blast, f"{fasta}_blast.tsv")

            print(f"Rodando BLASTp para {fasta}...")

            # --- Bloco de Comando Modificado ---
            comando = [
                "blastp",
                "-query", entrada,
                "-db", "pdb",
                "-remote",
                "-evalue", "1e-5",
                # Pega o valor do arquivo config
                "-max_target_seqs", str(config.blast_max_target_seqs), 
                "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
                "-out", saida
            ]
            # --- Fim do Bloco Modificado ---

            subprocess.run(comando)
            print(f"Resultado salvo em: {saida}")

        print("\nBLAST concluído!")
        
    except Exception as e:
        print(f"\nErro inesperado no BLAST: {e}")