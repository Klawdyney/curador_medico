import sqlite3

def preparar_banco():
    with sqlite3.connect('medical_insight.db') as conn:
        cursor = conn.cursor()
        # Lista de novas colunas para os Planos e Agendamento
        novas_colunas = [
            ("plano", "TEXT DEFAULT 'Básico'"),
            ("dia_envio", "TEXT DEFAULT 'Segunda'"),
            ("horario_envio", "TEXT DEFAULT '08:00'"),
            ("valor_assinatura", "REAL DEFAULT 89.90")
        ]
        
        for nome_col, tipo in novas_colunas:
            try:
                cursor.execute(f"ALTER TABLE clientes ADD COLUMN {nome_col} {tipo}")
                print(f"✅ Coluna {nome_col} adicionada com sucesso.")
            except sqlite3.OperationalError:
                print(f"⚠️ Coluna {nome_col} já existe (pulando...).")
    
    print("\n🚀 Banco de dados pronto para a Nova Fase!")

if __name__ == "__main__":
    preparar_banco()