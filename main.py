# main.py

import os
import sys
import config 

# Importa os módulos de funções
from pipeline_utils import filter_utils
from pipeline_utils import extract_utils
from pipeline_utils import blast_utils
from pipeline_utils import align_utils
from pipeline_utils import model_utils
from pipeline_utils import consensus_utils
from pipeline_utils import pdb_utils
from pipeline_utils import modeller_utils 

def encontrar_arquivo_input(diretorio, extensoes, tipo_arquivo):
    """
    Busca arquivos com extensões específicas na pasta 'input'.
    Se houver 1, retorna ele. Se houver >1, pede seleção.
    """
    if not os.path.exists(diretorio):
        return None

    arquivos = [f for f in os.listdir(diretorio) if f.lower().endswith(extensoes)]
    
    if not arquivos:
        print(f"[Aviso] Nenhum arquivo {tipo_arquivo} ({extensoes}) encontrado em 'input'.")
        return None
    
    if len(arquivos) == 1:
        print(f"[Input] Arquivo {tipo_arquivo} detectado: {arquivos[0]}")
        return os.path.join(diretorio, arquivos[0])
    
    print(f"\n[Input] Múltiplos arquivos {tipo_arquivo} encontrados. Qual utilizar?")
    for i, f in enumerate(arquivos, start=1):
        print(f"[{i}] {f}")
    
    while True:
        escolha = input("Escolha o número: ").strip()
        if escolha.isdigit():
            idx = int(escolha) - 1
            if 0 <= idx < len(arquivos):
                print(f"[Input] Selecionado: {arquivos[idx]}")
                return os.path.join(diretorio, arquivos[idx])
        print("Seleção inválida.")

def main():
    
    # --- 1. Configuração de Diretórios Base ---
    dir_pipeline = os.path.dirname(os.path.abspath(__file__)) 
    dir_input = os.path.join(dir_pipeline, "input")
    dir_results = os.path.join(dir_pipeline, "results")

    # --- 2. Definição dos Caminhos de Saída ---
    dir_f1 = os.path.join(dir_results, "Funcao1_Filtrar")
    dir_f2a = os.path.join(dir_results, "Funcao2a_Separar")
    dir_f2b = os.path.join(dir_results, "Funcao2b_FastasIndividuais")
    dir_f3 = os.path.join(dir_results, "Funcao3_Blastp")
    dir_f4 = os.path.join(dir_results, "Funcao4_AlinhamentoMultiplo")
    dir_f5 = os.path.join(dir_results, "Funcao5_Consensus")
    dir_f6 = os.path.join(dir_results, "Funcao6_PDB")
    dir_f7 = os.path.join(dir_results, "Funcao7_Modeller")

    pastas_output = [dir_f1, dir_f2a, dir_f2b, dir_f3, dir_f4, dir_f5, dir_f6, dir_f7]

    # --- 3. Criação de Diretórios ---
    os.makedirs(dir_input, exist_ok=True)
    os.makedirs(dir_results, exist_ok=True)
    for pasta in pastas_output:
        os.makedirs(pasta, exist_ok=True)
    
    print(f"Pipeline iniciado.")
    print(f"Pasta de Input: {dir_input}")
    print(f"Pasta de Resultados: {dir_results}")

    # --- 4. Detecção Automática de Inputs ---
    input_tsv = encontrar_arquivo_input(dir_input, (".tsv"), "InterPro TSV")
    input_fasta = encontrar_arquivo_input(dir_input, (".fasta"), "Sequências FASTA")

    while True:
        print("\nMENU PRINCIPAL--------------------------------------------------")
        print("1. Filtrar e Extrair Domínios      (Input: input; Output: Funcao1, Funcao2a)")
        print("2. Preparar FASTAs Individuais     (Input: Funcao1, Funcao2a, input; Output: Funcao2b)")
        print("3. Rodar BLASTp                    (Input: Funcao1, Funcao2a, Funcao2b, Funcao5; Output: Funcao3)")
        print("4. Alinhar domínios                (Input: Funcao1, Funcao2a; Output: Funcao4)")
        print("5. Gerar sequências consenso       (Input: Funcao4; Output: Funcao5)")
        print("6. Extrair e Baixar PDBs do BLAST  (Input: Funcao3; Output: Funcao6)")
        print("7. Rodar MODELLER                  (Input: Funcao2a, Funcao2b, Funcao6; Output: Funcao7)")
        print("------------------------------------------------------------------")

        opcao = input("Escolha uma opção (1-7): ").strip()
        
        # --- OPÇÃO 1: FILTRAR E EXTRAIR ---
        if opcao == '1':
            print(f"\n--- Iniciando Função 1: Filtrar e Extrair ---")
            
            if not input_tsv or not input_fasta:
                print("[ERRO] Arquivos de entrada não encontrados na pasta 'input'.")
                print("Certifique-se de colocar 1 arquivo .tsv e 1 arquivo .fasta lá.")
                
                retry = input("Deseja buscar novamente na pasta input? (s/n): ").strip().lower()
                if retry == 's':
                    input_tsv = encontrar_arquivo_input(dir_input, (".tsv"), "InterPro TSV")
                    input_fasta = encontrar_arquivo_input(dir_input, (".fasta", ".fa", ".fna"), "Sequências FASTA")
                    if not input_tsv or not input_fasta:
                        continue
                else:
                    continue

            print(f"\n[Parte 1] Filtrando proteínas...")
            df_filtrado, dominios_escolhidos, metodo_usado = filter_utils.filtrar_por_dominios_e_metodo(
                input_tsv, input_fasta, dir_f1
            )
            
            if df_filtrado is not None:
                print(f"\n[Parte 2a] Extraindo domínios...")
                extract_utils.extrair_outputs_fasta(
                    df_filtrado, dominios_escolhidos, metodo_usado, 
                    dir_f1, dir_f2a   
                )
                print("\nFunção 1 concluída! Arquivos salvos em 'Funcao1_Filtrar' e 'Funcao2a_Separar'.")
            else:
                print("\n[!] Filtragem falhou ou foi cancelada.")

        # --- OPÇÃO 2: PREPARAR FASTAS (SPLIT) ---
        elif opcao == '2':
            print(f"\n--- Iniciando Função 2: Preparar FASTAs Individuais ---")
            print("Onde estão os arquivos Multi-FASTA que você quer dividir?")
            print(f"[1] Pasta '{os.path.basename(dir_f2a)}' (Domínios extraídos)")
            print(f"[2] Pasta '{os.path.basename(dir_f1)}' (Proteínas filtradas)")
            print(f"[3] Pasta '{os.path.basename(dir_input)}' (Inserção Manual/Pasta Input)")
            
            escolha_fonte = input("Escolha (1, 2 ou 3): ").strip()
            
            source_dir = None
            if escolha_fonte == '1': source_dir = dir_f2a
            elif escolha_fonte == '2': source_dir = dir_f1
            elif escolha_fonte == '3': source_dir = dir_input
            else: 
                print("Opção inválida.")
                continue

            print(f"Lendo de: {source_dir}")
            print(f"Salvando individuais em: {dir_f2b}")
            
            model_utils.enviar_para_modelagem(source_dir, dir_f2b)

        # --- OPÇÃO 3: BLASTp ---
        elif opcao == '3':
            print(f"\n--- Iniciando Função 3: BLASTp ---")
            
            print("\nDe qual pasta você quer ler os arquivos FASTA?")
            print(f"[1] Da pasta '{os.path.basename(dir_f1)}' (Proteínas filtradas)")
            print(f"[2] Da pasta '{os.path.basename(dir_f2a)}' (Domínios separados)")
            print(f"[3] Da pasta '{os.path.basename(dir_f2b)}' (FASTAs individuais)")
            print(f"[4] Da pasta '{os.path.basename(dir_f5)}' (Sequências consenso)")
            
            escolha_fonte = input("Escolha (1, 2, 3 ou 4): ").strip()
            target_dir = dir_f3

            if escolha_fonte in ['1', '2']:
                # Modo Manual para pastas raiz (Pergunta 1 não existe, Pergunta 2 acontece)
                source_dir = dir_f1 if escolha_fonte == '1' else dir_f2a
                print(f"Lendo FASTAs de: {source_dir}")
                blast_utils.rodar_blast(source_dir, target_dir, automatico=False)

            elif escolha_fonte == '3':
                # Modo Automático para subpastas (Pergunta 1 acontece, Pergunta 2 suprimida)
                print(f"\nAnalisando subpastas em '{dir_f2b}'...")
                try:
                    subpastas = [d for d in os.listdir(dir_f2b) if os.path.isdir(os.path.join(dir_f2b, d))]
                    if not subpastas:
                        print("Nenhuma subpasta encontrada.")
                        continue
                    print("[0] TODAS as pastas abaixo")
                    for i, pasta in enumerate(subpastas, start=1): print(f"[{i}] {pasta}")
                    escolha = input("Escolha: ").strip()
                    pastas_proc = []
                    if escolha == '0': pastas_proc = subpastas
                    elif escolha.isdigit() and 1 <= int(escolha) <= len(subpastas):
                        pastas_proc.append(subpastas[int(escolha)-1])
                    for p in pastas_proc:
                        s_dir = os.path.join(dir_f2b, p)
                        t_dir = os.path.join(target_dir, p)
                        os.makedirs(t_dir, exist_ok=True)
                        blast_utils.rodar_blast(s_dir, t_dir, automatico=True)
                except Exception as e: print(f"Erro: {e}")

            elif escolha_fonte == '4':
                # Modo Automático para subpastas (Pergunta 1 acontece, Pergunta 2 suprimida)
                print(f"\nAnalisando subpastas em '{dir_f5}'...")
                try:
                    subpastas = [d for d in os.listdir(dir_f5) if os.path.isdir(os.path.join(dir_f5, d))]
                    if not subpastas:
                        print("Nenhuma subpasta encontrada.")
                        continue
                    print("[0] TODAS as pastas abaixo")
                    for i, pasta in enumerate(subpastas, start=1): print(f"[{i}] {pasta}")
                    escolha = input("Escolha: ").strip()
                    pastas_proc = []
                    if escolha == '0': pastas_proc = subpastas
                    elif escolha.isdigit() and 1 <= int(escolha) <= len(subpastas):
                        pastas_proc.append(subpastas[int(escolha)-1])
                    for p in pastas_proc:
                        s_dir = os.path.join(dir_f5, p)
                        # Padronização: nome da pasta de saída = nome da pasta de entrada (sem prefixo Consenso_)
                        t_dir = os.path.join(target_dir, p) 
                        os.makedirs(t_dir, exist_ok=True)
                        blast_utils.rodar_blast(s_dir, t_dir, automatico=True)
                except Exception as e: print(f"Erro: {e}")

        # --- OPÇÃO 4: ALINHAMENTO ---
        elif opcao == '4':
            print("\n--- Iniciando Função 4: Alinhamento ---")
            source_dir = ""
            while source_dir not in [dir_f1, dir_f2a]:
                escolha = input(f"Usar FASTAs da [1] '{os.path.basename(dir_f1)}' ou [2] '{os.path.basename(dir_f2a)}'? (1/2): ").strip()
                if escolha == '1': source_dir = dir_f1
                elif escolha == '2': source_dir = dir_f2a
                else: print("Seleção inválida.")
            print(f"Lendo FASTAs de: {source_dir}")
            print(f"Salvando em: {dir_f4}")
            align_utils.alinhar_dominios_clustalo_online(source_dir, dir_f4)
        
        # --- OPÇÃO 5: CONSENSO ---
        elif opcao == '5':
            print(f"\n--- Iniciando Função 5: Consenso ---")
            print(f"Lendo Alinhamentos de: {dir_f4}")
            print(f"Salvando em: {dir_f5}")
            consensus_utils.gerar_consensos_para_diretorio(dir_f4, dir_f5)
        
        # --- OPÇÃO 6: PDB ---
        elif opcao == '6':
            print(f"\n--- Iniciando Função 6: Extrair Códigos PDB ---")
            print(f"Lendo arquivos .tsv de: {dir_f3}")
            print(f"Salvando em: {dir_f6}")
            pdb_utils.extrair_pdb_codes(dir_f3, dir_f6)
       
        # --- OPÇÃO 7: MODELLER ---
        elif opcao == '7':
            print(f"\n--- Iniciando Função 7: Rodar MODELLER ---")
            # Adicionado dir_f5 na chamada para garantir a busca de consenso
            modeller_utils.run_modelling(dir_f2a, dir_f2b, dir_f5, dir_f6, dir_f7)
            
        else:
            print("\n[!] Opção inválida. Digite um número entre 1 e 7.")
            continue  

        continuar = input("\nDeseja realizar outra operação? (s/n): ").strip().lower()
        if continuar != 's': 
            print("\nEncerrando o programa... Obrigada por usar e até!! ;)") 
            break 

if __name__ == "__main__":
    main()