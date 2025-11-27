"""
Módulo para a Função 3: Alinhar domínios usando Clustal Omega (EBI).
"""

import os
import requests
import time
from Bio import AlignIO
from io import StringIO

def alinhar_dominios_clustalo_online(dir_leitura_fasta, dir_escrita_align):
    """
    Lê arquivos FASTA de 'Funcao1_Filtrar' ou 'Funcao2a_Separar' e salva os
    alinhamentos e árvores filogenéticas em subpastas
    dentro de 'Funcao4_AlinhamentoMultiplo'.
    """
    try:
        fasta_files = [f for f in os.listdir(dir_leitura_fasta) if f.endswith(".fasta")]
        if not fasta_files:
            print(f"\nNenhum arquivo FASTA encontrado em '{dir_leitura_fasta}'.")
            return None

        print("\nArquivos FASTA disponíveis para alinhamento:\n")
        for i, f in enumerate(fasta_files, start=1):
            print(f"[{i}] {f}")

        selecao = input(
            "\nDigite os números dos arquivos que deseja alinhar (ex: 1,3) ou 0 para todos: "
        ).strip()

        if selecao == "0":
            arquivos_escolhidos = fasta_files
        else:
            try:
                indices = [int(x.strip()) for x in selecao.split(",") if x.strip().isdigit()]
                arquivos_escolhidos = [fasta_files[i - 1] for i in indices if 1 <= i <= len(fasta_files)]
            except Exception:
                print("Seleção inválida.")
                return None

        if not arquivos_escolhidos:
            print("Nenhum arquivo selecionado. Encerrando.")
            return None

        email_usuario = input("Digite seu e-mail (ou 'n' para pular): ").strip()
        if email_usuario.lower() == 'n' or not email_usuario:
            email_usuario = "example@example.com"

        resultados = {}

        for arquivo_fasta in arquivos_escolhidos:
            print(f"\n--- Processando arquivo: {arquivo_fasta} ---")
            caminho_fasta = os.path.join(dir_leitura_fasta, arquivo_fasta)

            with open(caminho_fasta, 'r') as f:
                seq_data = f.read()

            url_run = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/"
            params = {
                'email': email_usuario,
                'stype': 'protein',
                'sequence': seq_data,
                'outfmt': 'clustal',
                'guidetreeout': 'true', 
            }

            print("Enviando para Clustal Omega...")
            response = requests.post(url_run, data=params)
            response.raise_for_status()
            job_id = response.text.strip()
            print(f"Job enviado! ID: {job_id}")

            url_status = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{job_id}"
            status = ""
            while status not in ["FINISHED", "ERROR"]:
                time.sleep(5)
                status_resp = requests.get(url_status)
                status_resp.raise_for_status()
                status = status_resp.text.strip()
                print(f"Status do job ({arquivo_fasta}): {status}")

            if status == "ERROR":
                print("Erro no job. Pulando arquivo.")
                continue

            # Organização de Pastas
            nome_base = os.path.splitext(arquivo_fasta)[0]
            pasta_saida_especifica = os.path.join(dir_escrita_align, nome_base)
            os.makedirs(pasta_saida_especifica, exist_ok=True)
            print(f"Salvando resultados em: {pasta_saida_especifica}")

            # Baixar e Salvar Alinhamento
            url_result = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/aln-clustal"
            result_resp = requests.get(url_result)
            result_resp.raise_for_status()
            aln_text = result_resp.text

            arquivo_alinhado = os.path.join(pasta_saida_especifica, f"{nome_base}_clustalo_alinhamento.clustal")
            with open(arquivo_alinhado, 'w') as f:
                f.write(aln_text)
            print(f"Alinhamento salvo em: {arquivo_alinhado}")

            alignment = AlignIO.read(StringIO(aln_text), "clustal")
            print(f"✓ {len(alignment)} sequências alinhadas ({alignment.get_alignment_length()} posições)")

            # Baixar e Salvar Árvore
            print("Baixando Phylogenetic Tree...")
            url_tree = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/tree"
            r = requests.get(url_tree)
            r.raise_for_status()
            
            arquivo_tree = os.path.join(pasta_saida_especifica, f"{nome_base}_tree.nwk")
            with open(arquivo_tree, 'w') as f:
                f.write(r.text)
            print(f"Árvore salva em: {arquivo_tree}")
            
            resultados[arquivo_fasta] = alignment

        print("\nTodos os alinhamentos concluídos com sucesso!")
        return resultados

    except Exception as e:
        print(f"\nErro no alinhamento Clustal Omega online: {e}")
        return None