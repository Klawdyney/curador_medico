from database import get_connection # Usa sua conexão que já funciona

def consertar_nuvem():
    print("🛠️ Adicionando colunas profissionais à nuvem...")
    conn = get_connection()
    cursor = conn.cursor()
    
    # As 3 colunas que o seu app.py está pedindo e a nuvem não tem
    colunas_faltantes = [
        ("status_assinatura", "TEXT DEFAULT 'ATIVO'"),
        ("dias_envio", "TEXT DEFAULT 'seg,qui'"),
        ("horario_envio", "TEXT DEFAULT '08:00'")
    ]
    
    for nome, tipo in colunas_faltantes:
        try:
            cursor.execute(f"ALTER TABLE clientes ADD COLUMN {nome} {tipo}")
            print(f"✅ Coluna {nome} adicionada!")
        except Exception:
            print(f"🟡 Coluna {nome} já existe. Ignorando.")
            
    conn.commit()
    conn.close()
    print("🏁 Nuvem sincronizada com o código profissional!")

if __name__ == "__main__":
    consertar_nuvem()