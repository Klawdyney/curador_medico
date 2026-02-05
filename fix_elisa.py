import sqlite3

# Conecta ao banco correto
conn = sqlite3.connect('medical_insight.db')
cursor = conn.cursor()

# Executa a atualização com os nomes de colunas e valores exatos que o robô precisa
cursor.execute("""
    UPDATE clientes 
    SET dia_envio = 'qua', 
        horario_envio = '14:00' 
    WHERE nome LIKE '%Elisa%'
""")

conn.commit()
print(f"✅ Sucesso! {cursor.rowcount} registro(s) da Elisa atualizado(s) para 'qua' às '14:00'.")
conn.close()