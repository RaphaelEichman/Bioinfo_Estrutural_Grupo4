"""
Módulo para a Função 1a: Filtrar por domínios e método.
"""

import os
import pandas as pd
from Bio import SeqIO
import config 

def filtrar_por_dominios_e_metodo(caminho_tsv, caminho_fasta, output_dir):
    """
    Filtra o TSV do InterPro e o FASTA de entrada.
    Salva os resultados em 'Funcao1_Filtrar'.
    """
    try:
        df = pd.read_csv(caminho_tsv, sep='\t', header=None)
        
        df[config.coluna_output] = df[config.coluna_output].astype(str)
        df[config.coluna_metodo] = df[config.coluna_metodo].astype(str)

        metodos_disponiveis = sorted(df[config.coluna_metodo].dropna().unique())
        print("\nMétodos de predição encontrados no arquivo:\n")
        for i, metodo in enumerate(metodos_disponiveis, start=1):
            print(f"{i}. {metodo}")

        while True:
            try:
                escolha = int(input("\nDigite o número do método que deseja utilizar: "))
                if 1 <= escolha <= len(metodos_disponiveis):
                    metodo_escolhido = metodos_disponiveis[escolha - 1]
                    break
                else:
                    print("Número inválido. Tente novamente.")
            except ValueError:
                print("Entrada inválida. Digite apenas o número correspondente ao método.")

        print(f"\nMétodo selecionado: {metodo_escolhido}")
        df_metodo = df[df[config.coluna_metodo].str.lower() == metodo_escolhido.lower()]

        outputs_disponiveis = sorted(df_metodo[config.coluna_output].dropna().unique())
        print("\nResultados disponíveis neste método:\n")
        for i, dom in enumerate(outputs_disponiveis, start=1):
            print(f"{i}. {dom}")

        outputs_input = input("\nDigite o(s) número(s) do(s) resultado(s) de seu interesse seguido de vírgulas (Ex. 1, 3, 4, 8): ").strip()
        outputs_de_interesse = []
        if all(item.strip().isdigit() for item in outputs_input.split(',')):
            indices = [int(i.strip()) for i in outputs_input.split(',') if i.strip().isdigit()]
            outputs_de_interesse = [outputs_disponiveis[i - 1] for i in indices if 1 <= i <= len(outputs_disponiveis)]

        if not outputs_de_interesse:
            print("Nenhum resultado foi informado. Encerrando.")
            return None, None, None

        print(f"\nResultados selecionados: {', '.join(outputs_de_interesse)}")

        busca = '|'.join(outputs_de_interesse)
        df_output = df_metodo[df_metodo[config.coluna_output].str.contains(busca, case=False, na=False)].copy()

        if df_output.empty:
            print(f"\nNenhuma proteína com os domínios {outputs_de_interesse} inferidos por {metodo_escolhido}.")
            return None, None, None

        ids_unicos = df_output[config.coluna_id].unique()
        resultados_sumarizados = []
        proteinas_todos_dominios = []

        for proteina_id in sorted(ids_unicos):
            df_proteina = df_output[df_output[config.coluna_id] == proteina_id]
            outputs_encontrados = [d for d in outputs_de_interesse
                                    if df_proteina[config.coluna_output].str.contains(d, case=False).any()]

            if outputs_encontrados:
                resultados_sumarizados.append({
                    'ID_Proteina': proteina_id,
                    'outputs_encontrados': ', '.join(outputs_encontrados)
                })

                if all(d in outputs_encontrados for d in outputs_de_interesse):
                    proteinas_todos_dominios.append(proteina_id)

        tsv_saida = os.path.join(output_dir, f"sumario_{metodo_escolhido}.tsv")
        df_final = pd.DataFrame(resultados_sumarizados)
        df_final.to_csv(tsv_saida, sep='\t', index=False)
        print(f"\n{len(df_final)} proteínas com os outputs de interesse ({metodo_escolhido}) foram encontradas.")
        print(f"Arquivo salvo em: '{tsv_saida}'")

        if proteinas_todos_dominios:
            arquivo_tres_dominios = os.path.join(output_dir, f"proteinas_{metodo_escolhido}_todos_outputs.txt")
            with open(arquivo_tres_dominios, 'w') as f:
                for p in proteinas_todos_dominios:
                    f.write(f"{p}\n")
            print(f"{len(proteinas_todos_dominios)} proteínas possuem todos os outputs informados.")
            print(f"Lista salva em: '{arquivo_tres_dominios}'")
        else:
            print("\nNenhuma proteína possui simultaneamente todos os domínios informados.")

        seq_dict = {record.id: record for record in SeqIO.parse(caminho_fasta, "fasta")}
        ids_filtrados = set(df_output[config.coluna_id].unique())

        arquivo_fasta_filtrado = os.path.join(output_dir, f"proteinas_filtradas_{metodo_escolhido}.fasta")
        with open(arquivo_fasta_filtrado, 'w') as out_fasta:
            for prot_id in ids_filtrados:
                if prot_id in seq_dict:
                    SeqIO.write(seq_dict[prot_id], out_fasta, "fasta")
        print(f"FASTA das proteínas filtradas salvo em: '{arquivo_fasta_filtrado}'")

        return df_output, outputs_de_interesse, metodo_escolhido

    except FileNotFoundError:
        print(f"\n[ERRO] Arquivo não encontrado. Verifique os caminhos:")
        print(f"TSV: {caminho_tsv}")
        print(f"FASTA: {caminho_fasta}")
        return None, None, None
    except Exception as e:
        print(f"\nOcorreu um erro: {e}")
        return None, None, None