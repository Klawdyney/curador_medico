import database_manager as db
from app import processar_medico_completo # Importa a função mestre validada
from datetime import datetime
import time
import logging
import sys

# Configuração de logs profissional
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout) # Garante que apareça no terminal da nuvem
    ]
)

def obter_dia_atual_sigla():
    """Converte o dia da semana para a sigla do banco (seg, ter, etc)."""
    dias = ['seg', 'ter', 'qua', 'qui', 'sex', 'sab', 'dom']
    return dias[datetime.now().weekday()]

def executar_rotina_agendada():
    """Lógica validada de verificação e envio."""
    dia_hoje = obter_dia_atual_sigla()
    # Pega a hora atual no formato exato do banco (ex: 14:00)
    hora_agora = datetime.now().strftime("%H:00")
    
    logging.info(f"🔎 Verificando agendamentos para {dia_hoje.upper()} às {hora_agora}...")
    
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
                
                # Chama a função MESTRA do app.py (PubMed -> Tradução -> Gemini -> PDF -> Envio)
                resultado = processar_medico_completo(medico)
                
                logging.info(f"🏁 Status Final: {resultado}")
                
                # Pausa técnica de 5s para não sobrecarregar APIs (Google/PubMed)
                time.sleep(5)
                
            except Exception as e:
                logging.error(f"❌ Erro ao processar {medico['nome']}: {e}")

    except Exception as e:
        logging.error(f"⚠️ Erro de conexão com o banco: {e}")

def iniciar_sentinela():
    """O Loop Infinito que mantém o robô vivo na nuvem."""
    print("\n" + "="*40)
    print("🤖 ROBÔ SENTINELA ATIVADO (Modo Contínuo)")
    print("👀 Monitorando relógio... Disparos apenas no minuto :00")
    print("="*40 + "\n")

    while True:
        agora = datetime.now()
        
        # O PULO DO GATO: Só trabalha se for o minuto 00 (Hora cheia)
        if agora.minute == 0:
            logging.info(f"⏰ HORA CHEIA DETECTADA ({agora.strftime('%H:%M')})! Acordando worker...")
            executar_rotina_agendada()
            
            # Dorme 65 segundos para garantir que saia do minuto 00 e não repita
            logging.info("💤 Ciclo concluído. Dormindo até a próxima hora...")
            time.sleep(65)
        
        else:
            # Se não for hora cheia, dorme o tempo que falta para o próximo minuto
            # Isso economiza CPU na nuvem e deixa o log limpo
            segundos_para_proximo_minuto = 60 - agora.second
            time.sleep(segundos_para_proximo_minuto)

if __name__ == "__main__":
    try:
        iniciar_sentinela()
    except KeyboardInterrupt:
        print("\n🛑 Sentinela desligado manualmente.")