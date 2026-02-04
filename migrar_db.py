import sqlite3

def migrar_banco_dados():
    # Conecta ao seu banco atual (verifique se o nome é exatamente clientes.db)
    conn = sqlite3.connect('clientes.db') 
    cursor = conn.cursor()

    # Novas colunas para o modelo SaaS (Fase 3)
    novas_colunas = [
        ("nome_clinica", "TEXT"),
        ("plano", "TEXT DEFAULT 'BASICO'"),
        ("status_assinatura", "TEXT DEFAULT 'ATIVO'"),
        ("dias_envio", "TEXT DEFAULT 'seg,qui'"),
        ("horario_envio", "TEXT DEFAULT '08:00'"),
        ("stripe_id", "TEXT"),
        ("ultimo_envio", "DATETIME")
    ]

    print("🚀 Iniciando migração do banco de dados...")

    for nome_coluna, tipo in novas_colunas:
        try:
            # Adiciona a coluna se ela não existir
            cursor.execute(f"ALTER TABLE clientes ADD COLUMN {nome_coluna} {tipo}")
            print(f"✅ Coluna '{nome_coluna}' adicionada.")
        except sqlite3.OperationalError:
            # Se a coluna já existir, ele apenas pula
            print(f"🟡 Coluna '{nome_coluna}' já existe. Mantendo dados.")

    conn.commit()
    conn.close()
    print("\n✨ Banco de dados atualizado com sucesso!")

if __name__ == "__main__":
    migrar_banco_dados()