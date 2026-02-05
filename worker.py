import database_manager as db
from app import processar_medico_completo # Importa a função mestre do seu app.py
from datetime import datetime
import time
import logging

# Configuração de logs para acompanhar o robô no terminal
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def obter_dia_atual_sigla():
    """Converte o dia da semana para a sigla do banco (seg, ter, etc)."""
    dias = ['seg', 'ter', 'qua', 'qui', 'sex', 'sab', 'dom']
    return dias[datetime.now().weekday()]

def executar_rotina_agendada():
    dia_hoje = obter_dia_atual_sigla()
    # Pega a hora atual no formato 10:00, 11:00, etc.
    hora_agora = datetime.now().strftime("%H:00")
    
    logging.info(f"🚀 Iniciando ciclo de curação para {dia_hoje.upper()} às {hora_agora}...")
    
    # 1. Busca todos os médicos ativos no Supabase
    medicos = db.buscar_todos_os_medicos_ativos()
    
    # 2. Filtra apenas os médicos agendados para este DIA e HORA
    # Isso encontrará o cadastro 'Claudinei' que você fez para as 10:00
    medicos_tarefa = [
        m for m in medicos 
        if m['dia_envio'] == dia_hoje and m['horario_envio'] == hora_agora
    ]
    
    if not medicos_tarefa:
        logging.info(f"📭 Nenhum médico agendado para o horário de {hora_agora}.")
        return

    logging.info(f"🩺 {len(medicos_tarefa)} médico(s) encontrado(s). Iniciando processamento...")

    for medico in medicos_tarefa:
        try:
            logging.info(f"🔄 Processando curadoria: Dr(a). {medico['nome']} ({medico['especialidade']})")
            
            # Chama a função que faz PubMed -> Gemini -> PDF -> Envio
            resultado = processar_medico_completo(medico)
            
            logging.info(f"Status: {resultado}")
            
            # Pausa de segurança para não sobrecarregar as APIs
            time.sleep(5)
            
        except Exception as e:
            logging.error(f"❌ Erro crítico ao processar {medico['nome']}: {e}")

if __name__ == "__main__":
    executar_rotina_agendada()