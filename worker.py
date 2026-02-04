import logging
from datetime import datetime, timedelta
from app import carregar_clientes_do_banco, processar_medico_completo

# Configura o log para aparecer no terminal do GitHub
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def executar_robot_curadoria():
    logging.info("🤖 Robo de Curadoria iniciado...")
    clientes = carregar_clientes_do_banco()
    
    # Ajuste de fuso: UTC para Brasil (Ipatinga/Timoteo)
    agora_brasil = datetime.now() - timedelta(hours=3)
    dias_map = {'mon':'seg','tue':'ter','wed':'qua','thu':'qui','fri':'sex','sat':'sab','sun':'dom'}
    dia_busca = dias_map.get(agora_brasil.strftime('%a').lower())
    hora_atual = agora_brasil.strftime('%H:00')

    logging.info(f"📅 Buscando: {dia_busca} as {hora_atual}")

    for id_c, user in clientes.items():
        if dia_busca in str(user['dias']).lower() and str(user['horario']) == hora_atual:
            try:
                logging.info(f"🚀 Iniciando motor para: {user['nome']}")
                resultado = processar_medico_completo(user)
                logging.info(f"🏁 Resultado: {resultado}")
            except Exception as e:
                logging.error(f"❌ Erro: {e}")

if __name__ == "__main__":
    executar_robot_curadoria()