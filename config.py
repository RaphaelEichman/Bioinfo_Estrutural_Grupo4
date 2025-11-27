# config.py
"""
Arquivo de configuração.
Armazena índices das colunas e parâmetros das ferramentas.
"""

# --- Colunas de interesse no tsv do InterProScan (sem header) ---
coluna_id = 0       # ID da protein
coluna_metodo = 3   # análise/DB
coluna_output = 5   # -o de acordo com a análise
coluna_start = 6    # start do -o
coluna_end = 7      # end do -o

# --- Configurações do BLAST ---
# Define o número máximo de sequências alvo (hits) 
# que o BLAST deve retornar.
blast_max_target_seqs = 10

# --- Configurações do MODELLER ---
# Define quantos modelos (PDBs) o MODELLER deve gerar.
modeller_ending_model = 5