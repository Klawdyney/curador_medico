import logging
from datetime import datetime, timedelta
from app import carregar_clientes_do_banco, processar_medico_completo

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def executar_robot_curadoria():
    logging.info("🤖 Robo de Curadoria iniciado...")
    clientes = carregar_clientes_do_banco()
    
    # Ajuste de fuso: Brasil (Ipatinga/Timoteo)
    agora_brasil = datetime.now() - timedelta(hours=3)
    dias_map = {'mon':'seg','tue':'ter','wed':'qua','thu':'qui','fri':'sex','sat':'sab','sun':'dom'}
    dia_busca = dias_map.get(agora_brasil.strftime('%a').lower())
    hora_atual = agora_brasil.strftime('%H:00')

    logging.info(f"📅 Relogio: {dia_busca} as {hora_atual} (Brasilia)")

    for id_c, user in clientes.items():
        # LOG DE SEGURANÇA: Mostra todos os dados do medico no terminal
        logging.info(f"📝 Dados lidos para {user.get('nome')}: {user}")

        # Tenta ler com os dois nomes possiveis (o antigo e o novo do banco)
        db_dia = str(user.get('dia_envio', user.get('dias', ''))).lower()
        db_hora = str(user.get('horario_envio', user.get('horario', ''))).strip()
        
        if dia_busca in db_dia and db_hora == hora_atual:
            try:
                logging.info(f"🚀 MATCH ENCONTRADO! Iniciando motor para: {user['nome']}")
                resultado = processar_medico_completo(user)
                logging.info(f"🏁 Resultado: {resultado}")
            except Exception as e:
                logging.error(f"❌ Erro no envio: {e}")
        else:
            logging.info(f"⏭️ Pulando {user['nome']}: Banco diz {db_dia}/{db_hora} mas agora e {dia_busca}/{hora_atual}")

if __name__ == "__main__":
    executar_robot_curadoria()