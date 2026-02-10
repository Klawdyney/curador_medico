import sqlite3

def criar_banco():
    conn = sqlite3.connect('medical_insight.db')
    cursor = conn.cursor()

    # Tabela de clientes (já está certa, mas mantemos aqui)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS clientes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nome TEXT NOT NULL,
            email TEXT NOT NULL UNIQUE,
            whatsapp TEXT,
            especialidade TEXT,
            clinica TEXT,
            keywords TEXT,
            plano TEXT,
            limite INTEGER DEFAULT 1,
            dia_envio TEXT,
            horario_envio TEXT
        )
    ''')

    # Tabela de histórico corrigida com PUBMED_ID
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS historico_envios (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            email_cliente TEXT,
            pubmed_id TEXT,
            titulo_artigo TEXT,
            data_envio TEXT,
            link_pubmed TEXT,
            FOREIGN KEY (email_cliente) REFERENCES clientes (email)
        )
    ''')

    conn.commit()
    conn.close()
    print("✅ Banco de dados atualizado com a coluna 'pubmed_id'!")

if __name__ == "__main__":
    criar_banco()