import logging
from datetime import datetime, timedelta # Adicionado timedelta para o fuso
# Importamos a 'Função Mestra' do app.py para garantir a mesma qualidade do PDF
from app import carregar_clientes_do_banco, processar_medico_completo

# Configuração de Logs para o terminal do GitHub Actions
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def executar_robot_curadoria():
    logging.info("🤖 Robô de Curadoria iniciado...")
    clientes = carregar_clientes_do_banco()
    
    # AJUSTE DE FUSO HORÁRIO: GitHub (UTC/Londres) para Brasil (UTC-3)
    # Isso garante que se você agendou 10:00, ele rode às 10:00 de Brasília
    agora_brasil = datetime.now() - timedelta(hours=3)
    
    dia_semana_hoje = agora_brasil.strftime('%a').lower() 
    # Mapeamento para bater com os nomes em português do seu banco
    dias_map = {'mon': 'seg', 'tue': 'ter', 'wed': 'qua', 'thu': 'qui', 'fri': 'sex', 'sat': 'sab', 'sun': 'dom'}
    dia_busca = dias_map.get(dia_semana_hoje)
    hora_atual = agora_brasil.strftime('%H:00')

    logging.info(f"📅 Verificando agendamentos para: {dia_busca} às {hora_atual} (Horário de Brasília)")

    for id_c, user in clientes.items():
        # LÓGICA DE FILTRO: Compara o agendamento do banco com o relógio de Brasília
        if dia_busca in str(user['dias']).lower() and user['horario'] == hora_atual:
            try:
                logging.info(f"🚀 Hora do show! Disparando motor para: {user['nome']}")
                
                # ACIONANDO A FUNÇÃO MESTRA: 
                # Aqui ele executa os 3 níveis, Gemini, PDF e Envio com a qualidade original
                resultado = processar_medico_completo(user)
                
                logging.info(f"Resultado: {resultado}")
            except Exception as e:
                logging.error(f"❌ Erro crítico no robô para {user['nome']}: {e}")

if __name__ == "__main__":
    executar_robot_curadoria()