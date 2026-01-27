import sqlite3

def configurar_historico():
    conexao = sqlite3.connect('medical_insight.db')
    cursor = conexao.cursor()
    # Cria a tabela de histórico vinculada ao ID do cliente
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS historico_envios (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            cliente_id TEXT,
            pubmed_id TEXT,
            data_envio DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (cliente_id) REFERENCES clientes (id)
        )
    ''')
    conexao.commit()
    conexao.close()
    print("Tabela de histórico criada com sucesso!")

if __name__ == "__main__":
    configurar_historico()