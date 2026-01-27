import sqlite3

def configurar_historico():
    try:
        conexao = sqlite3.connect('medical_insight.db')
        cursor = conexao.cursor()
        # Comando corrigido e completo
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS historico_envios (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                cliente_id TEXT,
                pubmed_id TEXT,
                data_envio DATETIME DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        conexao.commit()
        conexao.close()
        print("✅ SUCESSO: Tabela de histórico criada no banco de dados!")
    except Exception as e:
        print(f"❌ ERRO ao criar tabela: {e}")

if __name__ == "__main__":
    configurar_historico()