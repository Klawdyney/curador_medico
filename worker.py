import database_manager as db
from app import processar_medico_completo
from datetime import datetime, timedelta
import time
import logging
import sys
from concurrent.futures import ThreadPoolExecutor

# --- CONFIGURAÇÃO DE LOGS ---
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

def obter_dados_brasil():
    """
    Retorna o dia e a hora atuais ajustados para o Fuso Horário de Brasília (UTC-3).
    """
    # Pega a hora mundial (UTC) e subtrai 3 horas
    agora_brasil = datetime.utcnow() - timedelta(hours=3)
    
    dias_sigla = ['seg', 'ter', 'qua', 'qui', 'sex', 'sab', 'dom']
    dia_atual = dias_sigla[agora_brasil.weekday()]
    
    hora_atual = agora_brasil.strftime("%H:00")
    
    return dia_atual, hora_atual

def normalizar_dia_banco(dia_banco):
    """
    Traduz o que está no banco (ex: 'Quinta-feira') para a sigla que o sistema usa (ex: 'qui').
    """
    if not dia_banco: return ""
    
    dia_banco = dia_banco.lower()
    mapa = {
        'segunda-feira': 'seg', 'terça-feira': 'ter', 'quarta-feira': 'qua', 
        'quinta-feira': 'qui', 'sexta-feira': 'sex', 'sábado': 'sab', 'domingo': 'dom',
        'monday': 'seg', 'tuesday': 'ter', 'wednesday': 'qua', 'thursday': 'qui', 
        'friday': 'sex', 'saturday': 'sab', 'sunday': 'dom'
    }
    
    # Se já for sigla (ex: 'qui'), retorna ela mesma. Se for nome longo, traduz.
    return mapa.get(dia_banco, dia_banco[:3])

def tarefa_na_nuvem():
    logging.info("☁️ INICIANDO TAREFA (FUSO BRASIL -3h)...")
    
    # 1. Obtém hora e dia BRASIL
    dia_hoje, hora_agora = obter_dados_brasil()
    
    logging.info(f"📍 Relógio Brasil: {dia_hoje.upper()} às {hora_agora}")
    logging.info(f"🔎 Buscando agendamentos no banco...")
    
    try:
        medicos = db.buscar_todos_os_medicos_ativos()
        medicos_processar = []

        # 2. Filtragem Inteligente
        for m in medicos:
            # Traduz o dia do banco para sigla para comparar certo
            dia_medico = normalizar_dia_banco(m.get('dia_envio', ''))
            hora_medico = m.get('horario_envio', '')
            
            # Log de debug para vermos o que ele está lendo
            # logging.info(f"Checking: {m['nome']} | Dia: {dia_medico} vs {dia_hoje} | Hora: {hora_medico} vs {hora_agora}")
            
            if dia_medico == dia_hoje and hora_medico == hora_agora:
                medicos_processar.append(m)
        
        if not medicos_processar:
            logging.info(f"📭 Ninguém agendado para agora ({hora_agora}).")
            return
# 3. Processamento em Paralelo (Escala Profissional)
        logging.info(f"🩺 ENCONTRADO(S): {len(medicos_processar)} médicos para envio imediato.")

        try:
            # Processamos um por um com pausa para evitar o erro 429
            for medico in medicos_processar:
                processar_medico_completo(medico)
                time.sleep(2)  # Descanso vital de 2 segundos entre cada envio
                
            logging.info("✅ Ciclo de processamento sequencial concluído com sucesso.")
        except Exception as e:
            logging.error(f"❌ Erro durante o processamento: {e}")

    except Exception as e:
        logging.error(f"⚠️ Erro Geral no Worker: {e}")

if __name__ == "__main__":
    logging.info("🚀 Monitor de Escala Iniciado...")
    while True:
        tarefa_na_nuvem()
        # Agora ele acorda a cada 1 minuto para checar se mudou a hora
        logging.info("💤 Aguardando 20 minuto para próxima checagem...")
        time.sleep(1200)