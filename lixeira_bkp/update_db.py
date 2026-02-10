import sqlite3

def adicionar_coluna_whatsapp():
    conn = sqlite3.connect('medical_insight.db')
    cursor = conn.cursor()
    try:
        # Comando para adicionar a coluna nova
        cursor.execute('ALTER TABLE clientes ADD COLUMN whatsapp TEXT')
        conn.commit()
        print("✅ Sucesso: Coluna 'whatsapp' adicionada ao banco de dados!")
    except sqlite3.OperationalError:
        print("ℹ️ Aviso: A coluna 'whatsapp' já existe no banco.")
    finally:
        conn.close()

if __name__ == "__main__":
    adicionar_coluna_whatsapp()