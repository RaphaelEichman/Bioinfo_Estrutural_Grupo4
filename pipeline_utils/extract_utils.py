"""
Módulo para a Função 1b: Extrair sequências de todos os outputs.
"""

import os
from Bio import SeqIO
import config 

def extrair_outputs_fasta(df_output, outputs_de_interesse, metodo_escolhido, dir_leitura_fasta, dir_escrita_dominios):
    """
    Recebe os resultados da Função 1a, lê o FASTA filtrado de
    'Funcao1_Filtrar' e salva os domínios separados em 'Funcao2a_Separar'.
    """
    try:
        # Verifica se as variáveis de estado (da memória) são válidas
        if df_output is None or df_output.empty:
            print("\n[ERRO] Os dados de filtragem (df_output) estão vazios.")
            print("Isso pode acontecer se a etapa de filtragem não gerou resultados.")
            print("A execução da 'Função 1b' foi interrompida.")
            return

        print("\nLendo arquivo FASTA das proteínas filtradas...")
        
        # Constrói o caminho para o arquivo FASTA que a Função 1a deveria ter criado
        arquivo_fasta_filtrado = os.path.join(dir_leitura_fasta, f"proteinas_filtradas_{metodo_escolhido}.fasta")
        
        # Verifica se o arquivo FASTA filtrado realmente existe no disco
        if not os.path.exists(arquivo_fasta_filtrado):
            print(f"\n[ERRO] Arquivo FASTA filtrado não encontrado:")
            print(f"Caminho esperado: {arquivo_fasta_filtrado}")
            print("A 'Função 1a' pode ter falhado ao gerar este arquivo.")
            return
        
        # Se o arquivo existe, continua o processo
        seq_dict = {record.id: record for record in SeqIO.parse(arquivo_fasta_filtrado, "fasta")}

        for dominio in outputs_de_interesse:
            # --- LIMPEZA DE NOME DO ARQUIVO ---
            # Remove vírgulas e troca espaços por underscore
            output_name = dominio.replace(" ", "_").replace(",", "")
            # ----------------------------------

            # Salva o FASTA do domínio na pasta da Função 1b
            arquivo_fasta_output = os.path.join(dir_escrita_dominios, f"{output_name}_{metodo_escolhido}.fasta")
            
            # O parâmetro é 'case=False'
            df_dominio = df_output[df_output[config.coluna_output].str.contains(dominio, case=False, na=False)]

            if df_dominio.empty:
                print(f"\nNenhuma ocorrência do output '{dominio}' encontrada nos dados.")
                continue

            with open(arquivo_fasta_output, 'w') as fasta_out:
                for _, row in df_dominio.iterrows():
                    prot_id = row[config.coluna_id]
                    start = int(row[config.coluna_start])
                    end = int(row[config.coluna_end])

                    if prot_id not in seq_dict:
                        print(f"Aviso: ID {prot_id} encontrado no TSV mas não no FASTA filtrado.")
                        continue

                    seq_completa = seq_dict[prot_id].seq
                    sub_seq = seq_completa[start - 1:end]
                    
                    # Cria um header único: >ID_DominioLimpo_Inicio_Fim
                    header = f">{prot_id}_{output_name}_{start}_{end}"
                    fasta_out.write(f"{header}\n{sub_seq}\n")

            print(f"Sequências do output '{dominio}' salvas em: '{arquivo_fasta_output}'")
        print("\nExtração (Função 1b) concluída com sucesso.")

    except Exception as e:
        print(f"\nErro inesperado ao extrair as sequências dos domínios: {e}")