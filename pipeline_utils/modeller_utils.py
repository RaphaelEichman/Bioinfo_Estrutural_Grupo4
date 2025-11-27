# pipeline_utils/modeller_utils.py
"""
Módulo para a Função 7 (Modeller): Rodar o MODELLER.
(Versão Interativa: Seleção Automática via IDENTIDADE, Manual com Resolução)
"""

import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
import config
import sys
import contextlib
import json 
import re 

# Importações do MODELLER
try:
    from modeller import *
    from modeller.automodel import *
except ImportError:
    print("\n[ERRO] Modeller não encontrado. Verifique se está instalado e no seu PYTHONPATH.")
    
# Importa sys para usar float_info.max
import sys 

# -----------------------------------------------------------------
# CLASSE PERSONALIZADA DO AUTOMODEL
# -----------------------------------------------------------------
class MyAutoModel(AutoModel):
    def user_after_single_model(self):
        # Acessa o terminal original para dar um feedback limpo ao usuário
        sys.__stdout__.write(f"    -> Modelo gerado com sucesso.\n")
        sys.__stdout__.flush()

# -----------------------------------------------------------------
# FUNÇÃO HELPER: Encontrar a sequência alvo
# -----------------------------------------------------------------
def find_target_sequence(query_key, dir_f2b, dir_f2a):
    """Encontra um SeqRecord com um ID correspondente ao 'query_key'."""
    
    # 1. Busca em Funcao2b_FastasIndividuais
    for root, dirs, files in os.walk(dir_f2b):
        for file in files:
            caminho_fasta = os.path.join(root, file)
            file_name_without_ext = os.path.splitext(file)[0]

            if file_name_without_ext.startswith(query_key) and file.endswith(".fasta"):
                try:
                    record = SeqIO.read(caminho_fasta, "fasta")
                    if record.id == query_key: return record
                except ValueError:
                    try:
                        for record in SeqIO.parse(caminho_fasta, "fasta"):
                            if record.id == query_key: return record
                    except Exception: pass
                except Exception: pass

    # 2. Busca em Funcao2a_Separar (Domínios)
    try:
        for file in os.listdir(dir_f2a):
            if file.endswith(".fasta"):
                caminho_fasta = os.path.join(dir_f2a, file)
                try:
                    for record in SeqIO.parse(caminho_fasta, "fasta"):
                        if record.id == query_key: return record
                except Exception: pass
    except Exception: pass

    return None

# -----------------------------------------------------------------
# FUNÇÃO HELPER: Escrever o arquivo .ali
# -----------------------------------------------------------------
def write_sequence_to_ali(seq_record, ali_file_path, target_code="MtDH"):
    try:
        with open(ali_file_path, 'w') as f:
            f.write(f">P1;{target_code}\n")
            f.write(f"sequence:{target_code}:::::::0.00: 0.00\n")
            seq_str = str(seq_record.seq)
            f.write(f"{seq_str}*\n")
        return True
    except Exception as e:
        print(f"[ERRO] Falha ao criar arquivo .ali: {e}")
        return False

# -----------------------------------------------------------------
# FUNÇÃO HELPER: Extrair Resolução do PDB
# -----------------------------------------------------------------
def get_pdb_resolution(pdb_path):
    """
    Lê o arquivo PDB e tenta extrair a resolução (REMARK 2).
    Retorna uma string (ex: '1.50' ou 'N/A').
    """
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("REMARK   2 RESOLUTION."):
                    match = re.search(r"(\d+\.\d+)", line)
                    if match:
                        return match.group(1)
                    else:
                        return "N/A"
                if line.startswith("ATOM"):
                    break
    except Exception:
        return "Erro"
    
    return "N/A" 

# -----------------------------------------------------------------
# NOVO HELPER: Seleção Interativa de Molde
# -----------------------------------------------------------------
def selecionar_molde_interativo(run_dir):
    """
    Lê o JSON, verifica PDBs, extrai resolução e pede seleção.
    """
    json_path = os.path.join(run_dir, "blast_hits.json")
    pasta_moldes = os.path.join(run_dir, "Moldes")

    if not os.path.exists(json_path):
        print("  [ERRO] blast_hits.json não encontrado.")
        return None, None

    with open(json_path, 'r') as f:
        hits_data = json.load(f)
    
    available_hits = [] 
    
    for pdb_code, data in hits_data.items():
        pdb_full_path = os.path.join(pasta_moldes, f"{pdb_code}.pdb")
        if os.path.exists(pdb_full_path):
            
            resolution = get_pdb_resolution(pdb_full_path)
            
            available_hits.append({
                'code': pdb_code,
                'chain': data['chain'],
                'evalue': data.get('evalue', sys.float_info.max), 
                'bitscore': data.get('bitscore', 0),
                'pident': data.get('pident', 0),
                'resolution': resolution
            })

    if not available_hits:
        print("  [ERRO] Nenhum arquivo PDB válido encontrado na pasta Moldes.")
        return None, None

    # Prioridade: 1. Identidade (maior), 2. E-value (menor), 3. Bitscore (maior)
    available_hits.sort(key=lambda x: (-x['pident'], x['evalue'], -x['bitscore']))

    print(f"\n  --- Seleção de Molde para esta Proteína ---")
    print(f"  Moldes disponíveis: {len(available_hits)}")
    print("  Como deseja selecionar o molde?")
    print("  [1] Automático (Melhor Identidade)") 
    print("  [2] Manual (Ver lista completa)")
    
    while True:
        modo = input("  Escolha [1 ou 2]: ").strip()
        
        if modo == '1':
            melhor = available_hits[0]
            # --- LINHA ALTERADA AQUI ---
            print(f"  -> Selecionado Automaticamente: {melhor['code']} (Cadeia {melhor['chain']}) | Identidade: {melhor['pident']}% | Resolução: {melhor['resolution']} Å")
            # ---------------------------
            return melhor['code'], melhor['chain']
        
        elif modo == '2':
            print("\n  {:<5} {:<8} {:<8} {:<10} {:<12} {:<12} {:<12}".format(
                "ID", "PDB", "Cadeia", "Res.(Å)", "Ident.(%)", "Bitscore", "E-value"))
            print("  " + "-"*80)
            
            for i, hit in enumerate(available_hits, start=1):
                res_str = hit['resolution'] if hit['resolution'] else "N/A"
                print("  {:<5} {:<8} {:<8} {:<10} {:<12.2f} {:<12.1f} {:<12}".format(
                    f"[{i}]", 
                    hit['code'], 
                    hit['chain'], 
                    res_str, 
                    hit['pident'], 
                    hit['bitscore'], 
                    hit['evalue']
                ))
            
            print("  " + "-"*80)
            
            while True:
                sel = input(f"  Digite o número do molde (1-{len(available_hits)}): ").strip()
                try:
                    idx = int(sel) - 1
                    if 0 <= idx < len(available_hits):
                        escolhido = available_hits[idx]
                        print(f"  -> Selecionado Manualmente: {escolhido['code']} (Cadeia {escolhido['chain']})")
                        return escolhido['code'], escolhido['chain']
                    else:
                        print("  Número inválido.")
                except ValueError:
                    print("  Entrada inválida.")
        else:
            print("  Opção inválida. Digite 1 ou 2.")


# -----------------------------------------------------------------
# FUNÇÃO CORE DO MODELLER 
# -----------------------------------------------------------------
def _execute_modeller_core(run_dir, selected_code, selected_chain):
    """
    Executa a lógica interna do Modeller usando o molde JÁ SELECIONADO.
    """
    print(f"\n--- Iniciando Core do Modeller com {selected_code}:{selected_chain} ---")
    
    try:
        output_dir = "output" 
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)

        # Definição de variáveis locais para o Modeller
        X = selected_code
        X_CHAIN = selected_chain

        # --- ALINHAMENTO ---
        env_align = Environ()
        env_align.io.atom_files_directory = ['.', '../Moldes'] 
        
        aln_align = Alignment(env_align)
        
        # Carrega o molde selecionado
        mdl = Model(env_align, file=X, model_segment=('FIRST:'+X_CHAIN,'LAST:'+X_CHAIN))
        aln_align.append_model(mdl, align_codes=X+X_CHAIN, atom_files=X+'.pdb')
        
        # Carrega o alvo
        aln_align.append(file='../MtDH.ali', align_codes='MtDH') 
        
        # Alinha
        aln_align.align2d() 
        
        # Escreve arquivos de alinhamento
        aln_align.write(file=f'MtDH-{X}.ali', alignment_format='PIR')
        aln_align.write(file=f'MtDH-{X}.pap', alignment_format='PAP')

        # --- MODELAGEM ---
        env_model = Environ()
        env_model.io.atom_files_directory = ['.', '../Moldes'] 
        
        # Usamos MyAutoModel para feedback limpo no terminal
        a = MyAutoModel(env_model, alnfile=f'MtDH-{X}.ali',
                      knowns=X+X_CHAIN,
                      sequence='MtDH',
                      assess_methods=(assess.DOPE, assess.GA341))
        
        a.starting_model = 1
        a.ending_model = config.modeller_ending_model
        
        # Roda a modelagem
        a.make()
        
        return a.outputs

    except Exception as e:
        print(f"\n[ERRO NO CORE] {e}")
        raise e


# -----------------------------------------------------------------
# FUNÇÃO PRINCIPAL (ORQUESTRADOR)
# -----------------------------------------------------------------
def run_modelling(dir_f2a, dir_f2b, dir_f6, dir_f7):
    """
    Orquestra o processo de modelagem.
    """
    print("\n--- Iniciando Função 7: MODELLER ---")
    
    # --- 1. Selecionar Fonte (F6) ---
    try:
        pastas_base = [d for d in os.listdir(dir_f6) if os.path.isdir(os.path.join(dir_f6, d))]
    except Exception as e:
        print(f"Erro ao ler {dir_f6}: {e}")
        return

    if not pastas_base:
        print(f"Nenhuma pasta encontrada em {dir_f6}. Execute a Função 6 primeiro.")
        return

    print("Qual grupo de PDBs você deseja usar?")
    print("[0] TODAS as pastas abaixo") 
    for i, nome_pasta in enumerate(pastas_base, start=1):
        print(f"[{i}] {nome_pasta}")

    pastas_base_selecionadas = []
    try:
        escolha_base = input("Digite o número: ").strip()
        if escolha_base == "0":
            pastas_base_selecionadas = pastas_base
        elif 1 <= int(escolha_base) <= len(pastas_base):
            pastas_base_selecionadas.append(pastas_base[int(escolha_base) - 1])
        else:
            print("Seleção inválida.")
            return
    except ValueError:
        print("Entrada inválida.")
        return
    
    query_jobs_to_run = [] 
    
    for nome_pasta_base in pastas_base_selecionadas:
        dir_base_selecionada = os.path.join(dir_f6, nome_pasta_base)
        print(f"\nAnalisando pasta: {nome_pasta_base}")
        
        try:
            pastas_query_nomes = [d for d in os.listdir(dir_base_selecionada) if os.path.isdir(os.path.join(dir_base_selecionada, d))]
        except Exception as e:
            continue

        if not pastas_query_nomes:
            continue

        print("Qual proteína (query) você deseja modelar?")
        print("[0] TODAS as proteínas abaixo")
        for i, nome_pasta in enumerate(pastas_query_nomes, start=1):
            print(f"[{i}] {nome_pasta}")

        try:
            escolha_query = input(f"Digite o número para '{nome_pasta_base}': ").strip()
            
            if escolha_query == "0":
                for query_key in pastas_query_nomes:
                    template_source_path = os.path.join(dir_base_selecionada, query_key)
                    query_jobs_to_run.append((query_key, template_source_path, nome_pasta_base))
                
            elif 1 <= int(escolha_query) <= len(pastas_query_nomes):
                query_key = pastas_query_nomes[int(escolha_query) - 1]
                template_source_path = os.path.join(dir_base_selecionada, query_key)
                query_jobs_to_run.append((query_key, template_source_path, nome_pasta_base))
            else:
                print("Seleção inválida.")
        except ValueError:
            print("Entrada inválida.")
            
    if not query_jobs_to_run:
        print("\nNenhuma tarefa selecionada.")
        return
        
    print(f"\nIniciando {len(query_jobs_to_run)} tarefa(s)...")
    
    original_cwd = os.getcwd() 
    
    for query_key, template_source_path, nome_pasta_base in query_jobs_to_run:
        
        print(f"\n==================================================")
        print(f"PROCESSANDO: {query_key}")
        print(f"==================================================")
        
        # 2. Encontrar Alvo
        target_seq_record = find_target_sequence(query_key, dir_f2b, dir_f2a)
        if target_seq_record is None:
            print(f"[ERRO] Alvo não encontrado em F2b ou F2a. Pulando.")
            continue

        # 3. Preparar Estrutura de Pastas
        f7_run_base_dir = os.path.join(dir_f7, nome_pasta_base)
        os.makedirs(f7_run_base_dir, exist_ok=True)
        
        run_dir = os.path.join(f7_run_base_dir, query_key) 
        os.makedirs(run_dir, exist_ok=True)
        
        template_dest_path = os.path.join(run_dir, "Moldes")
        os.makedirs(template_dest_path, exist_ok=True)
        
        # Copiar JSON e PDBs
        json_source_path = os.path.join(template_source_path, "blast_hits.json")
        json_dest_path = os.path.join(run_dir, "blast_hits.json")
        
        if not os.path.exists(json_source_path):
            print(f"  [ERRO] 'blast_hits.json' não encontrado na origem. Pulando.")
            continue
        
        try:
            shutil.copy2(json_source_path, json_dest_path)
            for file in os.listdir(template_source_path):
                if file.endswith(".pdb"):
                    shutil.copy2(os.path.join(template_source_path, file),
                                os.path.join(template_dest_path, file))
        except Exception as e:
            print(f"  [ERRO] Falha ao copiar arquivos: {e}. Pulando.")
            continue
        
        # Criar .ali
        ali_file_path = os.path.join(run_dir, "MtDH.ali")
        if not write_sequence_to_ali(target_seq_record, ali_file_path, "MtDH"):
            continue

        # --- 4. SELEÇÃO DE MOLDE (Interativa ou Auto) ---
        selected_code, selected_chain = selecionar_molde_interativo(run_dir)
        
        if not selected_code:
            print("  [Abortado] Nenhum molde selecionado. Pulando esta proteína.")
            continue

        # --- 5. EXECUÇÃO DO MODELLER (Silenciada) ---
        log_file_path = os.path.join(run_dir, "saida.log")
        print(f"  Iniciando MODELLER... (Aguarde, log em: {os.path.basename(log_file_path)})")
        
        modeller_outputs = []
        
        try:
            with open(log_file_path, 'w') as log_file:
                # Redireciona TUDO para o arquivo
                with contextlib.redirect_stdout(log_file), contextlib.redirect_stderr(log_file):
                    os.chdir(run_dir)
                    modeller_outputs = _execute_modeller_core(run_dir, selected_code, selected_chain)
            
        except Exception as e:
            print(f"  [ERRO GERAL] Falha na execução do Modeller. Detalhes: {e}")
            os.chdir(original_cwd) 
            continue
        
        # --- 6. PÓS-PROCESSAMENTO (Pasta Selecionados) ---
        os.chdir(run_dir) 
        
        print("\n  --- Processando Melhores Resultados ---")
        dir_selecionados = os.path.join(run_dir, "Selecionados")
        os.makedirs(dir_selecionados, exist_ok=True)
        
        # A. Copiar Molde Selecionado
        caminho_molde_origem = os.path.join("Moldes", f"{selected_code}.pdb")
        caminho_molde_destino = os.path.join(dir_selecionados, f"Molde_{selected_code}.pdb")
        
        if os.path.exists(caminho_molde_origem):
            shutil.copy2(caminho_molde_origem, caminho_molde_destino)
            print(f"  -> Molde copiado: {os.path.basename(caminho_molde_destino)}")
        else:
            print(f"  -> [Aviso] Molde {selected_code}.pdb não encontrado para cópia.")

        # B. Copiar Melhor Modelo (Menor DOPE)
        try:
            output_folder = "output"
            if modeller_outputs:
                sucessos = [m for m in modeller_outputs if m.get('failure') is None]
                
                if sucessos:
                    melhor_modelo_data = min(sucessos, key=lambda x: x.get('DOPE score', 999999))
                    nome_arquivo_modelo = melhor_modelo_data['name'] 
                    valor_dope = melhor_modelo_data.get('DOPE score', 'N/A')
                    
                    caminho_modelo_origem = os.path.join(output_folder, nome_arquivo_modelo)
                    
                    if os.path.exists(caminho_modelo_origem):
                        extensao = os.path.splitext(nome_arquivo_modelo)[1]
                        nome_destino = f"Melhor_Modelo_DOPE_{valor_dope:.2f}{extensao}"
                        shutil.copy2(caminho_modelo_origem, os.path.join(dir_selecionados, nome_destino))
                        print(f"  -> Melhor modelo copiado: {nome_destino} (DOPE: {valor_dope})")
                    else:
                        print(f"  -> [Erro] Arquivo {nome_arquivo_modelo} não encontrado em {output_folder}.")
                else:
                    print("  -> Nenhum modelo gerado com sucesso.")
            else:
                print("  -> Lista de outputs do Modeller vazia.")
                
        except Exception as e:
            print(f"  -> Erro ao copiar melhores resultados: {e}")

        os.chdir(original_cwd)
            
    print(f"\n--- Função 7 (MODELLER) concluída ---")