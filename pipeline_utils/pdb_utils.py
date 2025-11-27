# pipeline_utils/pdb_utils.py
"""
Módulo para a Função 5: Extrair códigos PDB, salvar scores do BLAST 
e baixar os arquivos .pdb, agrupados por proteína (query).
"""

import os
import pandas as pd
import requests 
import time     
import json 

# -----------------------------------------------------------------
# FUNÇÃO AUXILIAR: Baixar PDBs
# (Sem alterações)
# -----------------------------------------------------------------
def baixar_pdb_files(codigos_pdb_set, pasta_saida_especifica):
    """
    Baixa uma lista/set de códigos PDB para uma pasta de saída específica.
    """
    try:
        os.makedirs(pasta_saida_especifica, exist_ok=True)
        
        if not codigos_pdb_set:
            print("  -> Nenhum código PDB para baixar neste grupo.")
            return

        print(f"  -> Baixando {len(codigos_pdb_set)} PDBs para '{os.path.basename(pasta_saida_especifica)}'...")
        
        for code in sorted(list(codigos_pdb_set)):
            code = code.strip().upper()
            if not code:
                continue
            
            url = f"https://files.rcsb.org/download/{code}.pdb"
            caminho_pdb_out = os.path.join(pasta_saida_especifica, f"{code}.pdb")

            if os.path.exists(caminho_pdb_out):
                print(f"    -> {code}.pdb já existe. Pulando.")
                continue

            try:
                r = requests.get(url, timeout=10)
                r.raise_for_status() 
                
                with open(caminho_pdb_out, 'w') as f_pdb:
                    f_pdb.write(r.text)
                print(f"    -> {code}.pdb baixado com sucesso.")

                time.sleep(0.5) 

            except requests.exceptions.RequestException as e:
                print(f"    -> Falha ao baixar {code}.pdb (Erro: {e})")

        print("  -> Download de PDBs para este grupo concluído.")

    except Exception as e:
        print(f"  -> Erro inesperado durante o download de PDBs: {e}")

# -----------------------------------------------------------------
# FUNÇÃO PRINCIPAL (MODIFICADA)
# -----------------------------------------------------------------
def extrair_pdb_codes(dir_leitura_blast, dir_escrita_pdb):
    """
    Varre 'dir_leitura_blast', agrupa hits por proteína (Coluna A),
    salva um 'blast_hits.json' com os scores E A CADEIA, e baixa os PDBs.
    """
    
    arquivos_tsv_encontrados = []
    
    print(f"Buscando arquivos .tsv em: {dir_leitura_blast}...")
    for root, dirs, files in os.walk(dir_leitura_blast):
        for file in files:
            if file.endswith(".tsv"):
                arquivos_tsv_encontrados.append(os.path.join(root, file))
    
    if not arquivos_tsv_encontrados:
        print(f"Nenhum arquivo .tsv encontrado em '{dir_leitura_blast}'.")
        print("Execute a Função 2 (BLASTp) primeiro.")
        return

    print(f"Encontrados {len(arquivos_tsv_encontrados)} arquivos .tsv para processar...\n")
    
    total_proteinas_processadas = 0

    for caminho_tsv in arquivos_tsv_encontrados:
        arquivo_base = os.path.basename(caminho_tsv)
        nome_base_pasta_tsv = os.path.splitext(arquivo_base)[0]
        print(f"Processando arquivo: {arquivo_base}")

        pasta_saida_base_tsv = os.path.join(dir_escrita_pdb, nome_base_pasta_tsv)
        os.makedirs(pasta_saida_base_tsv, exist_ok=True)
        
        try:
            df = pd.read_csv(caminho_tsv, sep='\t', header=None, on_bad_lines='skip')
            
            if df.empty:
                print(f"  -> Aviso: Arquivo {arquivo_base} não contém dados. Pulando.")
                continue 

            # Colunas: 0=Query, 1=Hit, 2=P-ident, 4=E-value, 5=Bitscore
            if 5 not in df.columns:
                print(f"  -> Aviso: Arquivo {arquivo_base} não parece ter 6+ colunas. Pulando.")
                continue
            
            grupos_de_proteinas = df.groupby(0)
            
            print(f"  -> Encontradas {len(grupos_de_proteinas)} proteínas (queries) neste arquivo.")
            total_proteinas_processadas += len(grupos_de_proteinas)

            for query_id, group_df in grupos_de_proteinas:
                
                nome_pasta_proteina = "".join(
                    c for c in query_id if c.isalnum() or c in ('_', '-')
                ).rstrip()
                if len(nome_pasta_proteina) > 100: 
                    nome_pasta_proteina = nome_pasta_proteina[:100]

                print(f"\n    Processando Query: {query_id}")

                pasta_saida_proteina = os.path.join(pasta_saida_base_tsv, nome_pasta_proteina)
                os.makedirs(pasta_saida_proteina, exist_ok=True)
                
                codigos_neste_grupo = set()
                hits_data = {} 

                for index, row in group_df.iterrows():
                    hit_string = str(row[1])
                    parts = hit_string.split('|')
                    
                    # --- LÓGICA DE EXTRAÇÃO ATUALIZADA ---
                    # Formato esperado: pdb|CODIGO|CADEIA
                    if len(parts) >= 3 and parts[0] == 'pdb': 
                        pdb_code = parts[1]
                        chain = parts[2] # <-- Salva a cadeia
                        codigos_neste_grupo.add(pdb_code)
                        
                        if (pdb_code not in hits_data) or (row[4] < hits_data[pdb_code]["evalue"]):
                            hits_data[pdb_code] = {
                                "chain": chain, # <-- Salva a cadeia no JSON
                                "evalue": row[4],
                                "bitscore": row[5],
                                "pident": row[2]
                            }
                    # --- FIM DA LÓGICA ATUALIZADA ---
                
                if codigos_neste_grupo:
                    json_path = os.path.join(pasta_saida_proteina, 'blast_hits.json')
                    try:
                        with open(json_path, 'w') as f:
                            json.dump(hits_data, f, indent=4)
                        print(f"    -> 'blast_hits.json' (com cadeias) salvo em '{os.path.basename(pasta_saida_proteina)}'")
                    except Exception as e:
                        print(f"    -> Erro ao salvar JSON: {e}")
                    
                    baixar_pdb_files(codigos_neste_grupo, pasta_saida_proteina)
                else:
                    print(f"    -> Nenhum código PDB no formato pdb|CODIGO|CADEIA encontrado.")

        except pd.errors.EmptyDataError:
            print(f"  -> Aviso: Arquivo {arquivo_base} está vazio. Pulando.")
        except Exception as e:
            print(f"  -> Erro ao processar {arquivo_base}: {e}")

    print(f"\nExtração e download de PDB concluídos. {total_proteinas_processadas} queries processadas.")
    print(f"Resultados salvos em subpastas dentro de: {dir_escrita_pdb}\n")