# pipeline_utils/consensus_utils.py
"""
Módulo para a Função 4: Gerar sequências consenso de alinhamentos.
"""

from Bio import AlignIO
from collections import Counter
import random
import os

# ---------------------------------------------------------------
# Dicionário de classes de aminoácidos
# (Sem alterações)
# ---------------------------------------------------------------
residue_classes = {
    "A": "hidrofóbico", "V": "hidrofóbico", "L": "hidrofóbico", "I": "hidrofóbico", "M": "hidrofóbico",
    "F": "aromático", "Y": "aromático", "W": "aromático",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar",
    "K": "carregado positivo", "R": "carregado positivo", "H": "carregado positivo",
    "D": "carregado negativo", "E": "carregado negativo",
    "C": "especial (enxofre)", "G": "especial (flexível)", "P": "especial (cíclico)"
}

# ---------------------------------------------------------------
# Função para processar um único arquivo (lógica interna)
# (Sem alterações - esta função já salva onde é mandada)
# ---------------------------------------------------------------
def gerar_consenso_e_relatorio(arquivo_alinhamento, formato="clustal", limite_gaps=0.7, pasta_saida="consensus_results"):
    """
    Processa um único arquivo de alinhamento e salva seu consenso
    e relatório na 'pasta_saida' especificada.
    """
    alinhamento = AlignIO.read(arquivo_alinhamento, formato)
    n_seq = len(alinhamento)
    tamanho = alinhamento.get_alignment_length()
    
    nome_base = os.path.splitext(os.path.basename(arquivo_alinhamento))[0]
    
    consenso = []
    relatorio = []

    for pos in range(tamanho):
        coluna = [seq[pos] for seq in alinhamento]
        contagem = Counter(coluna)
        total = sum(contagem.values())
        
        sem_gaps = {aa: freq for aa, freq in contagem.items() if aa != "-"}
        gaps = contagem.get("-", 0)
        freq_gaps = gaps / total
        
        alerta = ""
        if not sem_gaps:
            consenso.append("-")
            alerta = "Posição totalmente com gaps"
            freq_consenso = 0
            aa_consenso = "-"
        else:
            max_freq = max(sem_gaps.values())
            top_residuos = [aa for aa, f in sem_gaps.items() if f == max_freq]
            aa_consenso = random.choice(top_residuos) if len(top_residuos) > 1 else top_residuos[0]
            freq_consenso = (max_freq / total) * 100
            
            if len(top_residuos) > 1:
                alerta = f"Empate entre {', '.join(top_residuos)}"
            elif freq_consenso < 50:
                alerta = f"{aa_consenso} presente em apenas {freq_consenso:.1f}% das sequências"
            elif freq_gaps > limite_gaps:
                alerta = f"Região com {freq_gaps*100:.1f}% de gaps"
        
            consenso.append(aa_consenso)
        
        classe = residue_classes.get(aa_consenso, "desconhecida")
        relatorio.append([
            pos + 1, aa_consenso, f"{freq_consenso:.1f}", n_seq, classe, alerta
        ])
    
    os.makedirs(pasta_saida, exist_ok=True)
    
    consenso_seq = "".join(consenso)
    
    fasta_saida = os.path.join(pasta_saida, f"{nome_base}_consensus.fasta")
    relatorio_saida = os.path.join(pasta_saida, f"{nome_base}_report.tsv")
    
    with open(fasta_saida, "w") as f:
        f.write(f">Consensus_{nome_base}\n")
        for i in range(0, len(consenso_seq), 60):
            f.write(consenso_seq[i:i+60] + "\n")
    
    with open(relatorio_saida, "w") as f:
        f.write("Posição\tResíduo_Consenso\tFrequência_Consenso(%)\tNº_Sequências\tClasse_Residuo\tAlerta\n")
        for linha in relatorio:
            f.write("\t".join(map(str, linha)) + "\n")
    
    print(f"  -> Salvo em: {pasta_saida}")
    
    return f"{nome_base}_consensus", consenso_seq

# ---------------------------------------------------------------
# Função principal (MODIFICADA)
# ---------------------------------------------------------------
def gerar_consensos_para_diretorio(dir_leitura_align, dir_escrita_consensus, limite_gaps=0.7):
    """
    Busca (recursivamente) por arquivos de alinhamento em 'dir_leitura_align',
    gera consenso e salva os resultados em subpastas espelhadas
    dentro de 'dir_escrita_consensus'.
    """
    # Esta é a pasta RAIZ (Funcao4_Consensus)
    pasta_saida_raiz = dir_escrita_consensus 

    # --- LÓGICA DE BUSCA MODIFICADA ---
    # Lista para guardar tuplas: (caminho_completo_input, pasta_saida_especifica_output)
    arquivos_a_processar = [] 
    
    print(f"Buscando arquivos de alinhamento em: {dir_leitura_align}...")
    
    # os.walk() irá percorrer o dir_leitura_align e todas as suas subpastas
    for root, dirs, files in os.walk(dir_leitura_align):
        
        # Encontra o caminho relativo da subpasta
        # Ex: "Diguanylate_cyclase" ou "." (se estiver na raiz de F3)
        rel_path = os.path.relpath(root, dir_leitura_align)
        
        # Define a pasta de saída espelhada
        # Ex: ".../Funcao4_Consensus/Diguanylate_cyclase"
        if rel_path == ".":
            pasta_saida_especifica = pasta_saida_raiz
        else:
            pasta_saida_especifica = os.path.join(pasta_saida_raiz, rel_path)

        # Encontra os arquivos de alinhamento dentro desta subpasta
        for file in files:
            if file.endswith((".clustal", ".aln", ".fasta")):
                caminho_completo_input = os.path.join(root, file)
                
                # Armazena o par de (input, output)
                arquivos_a_processar.append((caminho_completo_input, pasta_saida_especifica))
    
    if not arquivos_a_processar:
        print(f"Nenhum arquivo de alinhamento (.clustal, .aln, .fasta) encontrado em '{dir_leitura_align}' ou suas subpastas.")
        print("Execute a Função 3 primeiro.")
        return

    print(f"Encontrados {len(arquivos_a_processar)} arquivos de alinhamento para processar...\n")

    todas_consensos = [] 

    # Itera sobre a lista de pares (input, output)
    for caminho_completo, pasta_saida_especifica in arquivos_a_processar:
        
        # Garante que a subpasta de saída (ex: .../F4/Diguanylate_cyclase) exista
        os.makedirs(pasta_saida_especifica, exist_ok=True) 
        
        arquivo_base = os.path.basename(caminho_completo)
        formato = "clustal" if arquivo_base.endswith((".clustal", ".aln")) else "fasta"
        
        print(f"Processando: {arquivo_base}")
        try:
            # Passa o caminho completo do arquivo
            # E a pasta de SAÍDA específica (a subpasta)
            nome_base, consenso_seq = gerar_consenso_e_relatorio(
                caminho_completo, 
                formato=formato, 
                limite_gaps=limite_gaps, 
                pasta_saida=pasta_saida_especifica 
            )
            todas_consensos.append((nome_base, consenso_seq))
        except Exception as e:
            print(f"Erro ao processar {arquivo_base}: {e}")

    # --- Salva o arquivo 'todas_consensus.fasta' na pasta RAIZ ---
    if todas_consensos:
        # 'pasta_saida_raiz' é a 'Funcao4_Consensus'
        fasta_geral = os.path.join(pasta_saida_raiz, "todas_consensus.fasta")
        with open(fasta_geral, "w") as f:
            for nome_base, seq in todas_consensos:
                f.write(f">Consensus_{nome_base}\n") 
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")
        print(f"\n Arquivo geral criado: {fasta_geral}")

    print("\nProcessamento de consenso concluído!")
    print(f"Resultados individuais salvos em subpastas de: {pasta_saida_raiz}\n")