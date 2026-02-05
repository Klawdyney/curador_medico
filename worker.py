import database_manager as db
from app import processar_medico_completo
from datetime import datetime
import time
import logging
import sys

# --- CONFIGURAÇÃO DE LOGS ---
# Ajustado para exibir logs no painel do GitHub Actions
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

def obter_dia_atual_sigla():
    """Converte o dia da semana para a sigla do banco (seg, ter, etc)."""
    dias = ['seg', 'ter', 'qua', 'qui', 'sex', 'sab', 'dom']
    return dias[datetime.now().weekday()]

def tarefa_na_nuvem():
    """
    MODO NUVEM: Esta função roda, verifica os agendamentos da hora cheia,
    processa os envios e FINALIZA. 
    Não usamos 'while True' aqui para não travar o servidor do GitHub.
    """
    print("☁️ INICIANDO TAREFA AGENDADA NA NUVEM...")
    
    dia_hoje = obter_dia_atual_sigla()
    
    # Pega a hora cheia atual (ex: Se rodar às 14:05, pega "14:00")
    # Isso garante o sincronismo com o agendamento do GitHub
    hora_agora = datetime.now().strftime("%H:00")
    
    logging.info(f"🔎 Verificando envios para {dia_hoje.upper()} às {hora_agora}...")
    
    try:
        # 1. Busca todos os médicos ativos no Supabase
        medicos = db.buscar_todos_os_medicos_ativos()
        
        # 2. Filtra apenas os médicos agendados para este DIA e HORA
        medicos_tarefa = [
            m for m in medicos 
            if m['dia_envio'] == dia_hoje and m['horario_envio'] == hora_agora
        ]
        
        if not medicos_tarefa:
            logging.info(f"📭 Nenhum envio programado para agora ({hora_agora}).")
            return

        logging.info(f"🩺 ENCONTRADO(S): {len(medicos_tarefa)} médico(s). Iniciando processamento...")

        for medico in medicos_tarefa:
            try:
                logging.info(f"🚀 Iniciando curadoria: Dr(a). {medico['nome']} ({medico['especialidade']})")
                
                # Chama a função MESTRA validada do app.py
                resultado = processar_medico_completo(medico)
                
                logging.info(f"🏁 Status Final: {resultado}")
                
                # Pausa técnica de 5s para não sobrecarregar APIs
                time.sleep(5)
                
            except Exception as e:
                logging.error(f"❌ Erro ao processar {medico['nome']}: {e}")

    except Exception as e:
        logging.error(f"⚠️ Erro de conexão com o banco ou processamento geral: {e}")

if __name__ == "__main__":
    # Removemos o 'iniciar_sentinela' e o 'while True'.
    # Agora ele roda uma vez e encerra, perfeito para automação.
    tarefa_na_nuvem()