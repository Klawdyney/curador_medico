import logging
from datetime import datetime, timedelta
from app import carregar_clientes_do_banco, processar_medico_completo

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def executar_robot_curadoria():
    logging.info("🤖 Robo de Curadoria iniciado...")
    clientes = carregar_clientes_do_banco()
    
    # Ajuste de fuso: Brasil
    agora_brasil = datetime.now() - timedelta(hours=3)
    dias_map = {'mon':'seg','tue':'ter','wed':'qua','thu':'qui','fri':'sex','sat':'sab','sun':'dom'}
    dia_busca = dias_map.get(agora_brasil.strftime('%a').lower())
    hora_atual = agora_brasil.strftime('%H:00')

    logging.info(f"📅 Relógio: {dia_busca} às {hora_atual} (Brasília)")

    for id_c, user in clientes.items():
        # AQUI ESTÁ A CORREÇÃO: Usando os nomes que aparecem na sua imagem do SQL
        db_dia = str(user.get('dia_envio', '')).lower()
        db_hora = str(user.get('horario_envio', '')).strip()
        
        logging.info(f"🔍 Medico: {user['nome']} | Banco diz: {db_dia} às {db_hora}")

        if dia_busca in db_dia and db_hora == hora_atual:
            try:
                logging.info(f"🚀 CONDIÇÃO ACEITA! Enviando para: {user['nome']}")
                resultado = processar_medico_completo(user)
                logging.info(f"✅ Resultado: {resultado}")
            except Exception as e:
                logging.error(f"❌ Erro no envio: {e}")