import sqlite3

def atualizar_horario_teste():
    # Caminho corrigido para o seu banco real
    conn = sqlite3.connect('medical_insight.db')
    cursor = conn.cursor()

    # Ajustado para 12:00 para o próximo teste
    novo_horario = "12:00" 
    novo_dia = "qua"

    try:
        cursor.execute('''
            UPDATE clientes 
            SET horario_envio = ?, dia_envio = ? 
            WHERE nome LIKE '%Claudinei%'
        ''', (novo_horario, novo_dia))

        if cursor.rowcount > 0:
            conn.commit()
            print(f"✅ Sucesso! Horário do Dr. Claudinei atualizado para {novo_horario}.")
        else:
            print("⚠️ Médico não encontrado. Verifique se o nome está correto no banco.")
            
    except Exception as e:
        print(f"❌ Erro ao acessar a tabela: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    atualizar_horario_teste()