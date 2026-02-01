import streamlit as st
import sqlite3

# Configuração da Página
st.set_page_config(page_title="Medical In-Sight | Portal do Assinante", page_icon="🩺", layout="wide")

# Estilo CSS Personalizado
st.markdown("""
    <style>
    .main { background-color: #f8f9fa; }
    .stButton>button { width: 100%; border-radius: 5px; height: 3em; background-color: #003366; color: white; }
    .stExpander { background-color: white; border-radius: 10px; box-shadow: 0px 2px 10px rgba(0,0,0,0.05); border: 1px solid #e0e0e0; margin-bottom: 10px; }
    </style>
    """, unsafe_allow_html=True)

# --- FUNÇÕES DE BANCO DE DADOS ---

def buscar_dados_medico(email):
    try:
        conn = sqlite3.connect('medical_insight.db')
        cursor = conn.cursor()
        # Buscamos o que o cadastro básico oferece
        cursor.execute("SELECT nome, especialidade, keywords, plano FROM clientes WHERE email = ?", (email,))
        resultado = cursor.fetchone()
        conn.close()
        return resultado
    except:
        return None

def buscar_historico(email):
    """Nova função para buscar os artigos já enviados"""
    try:
        conn = sqlite3.connect('medical_insight.db')
        cursor = conn.cursor()
        cursor.execute("""
            SELECT titulo_artigo, data_envio, link_pubmed 
            FROM historico_envios 
            WHERE email_cliente = ? 
            ORDER BY id DESC
        """, (email,))
        rows = cursor.fetchall()
        conn.close()
        return rows
    except Exception as e:
        return []

def exibir_faq_completo():
    st.write("### 🛡️ Central de Transparência e Ajuda")
    st.info("Nossa operação é baseada em evidência científica e ética de dados. Veja como trabalhamos:")

    with st.expander("🔍 1. VERACIDADE E QUALIDADE DAS INFORMAÇÕES"):
        st.write("""
        **Como o Medical In-Sight garante que as informações são verídicas?**
        Nosso sistema conecta-se diretamente à API do PubMed/MEDLINE. Não utilizamos blogs ou sites de notícias; apenas artigos científicos com DOI rastreável. Cada boletim inclui o link direto para a fonte original.

        **Existe risco de 'alucinação' da IA nos relatórios?**
        Implementamos o protocolo 'Grounded AI'. Nosso motor é restrito aos fatos do artigo técnico. Utilizamos filtros de rigor metodológico para priorizar estudos de alto impacto (Ensaios Clínicos e Metanálises).
        """)

    with st.expander("🎯 2. OBJETIVO E FINALIDADE DO SERVIÇO"):
        st.write("""
        **O Medical In-Sight substitui a minha leitura técnica?**
        De forma alguma. Nossa finalidade é a curadoria e otimização de tempo. Preparamos o terreno para que você dedique seu tempo apenas ao que realmente impacta sua conduta.

        **O que acontece se não houver novidades na minha área em uma semana?**
        Prezamos pela qualidade. Se não houver publicações de alto impacto recentes, nosso Radar Positivo seleciona um 'Estudo Clássico' (Landmark Trial) fundamental para sua especialidade.
        """)

    with st.expander("🔐 3. SEGURANÇA E PRIVACIDADE (LGPD)"):
        st.write("""
        **Meus dados estão seguros?**
        Sim. Operamos sob os pilares da LGPD. Seus dados e keywords são criptografados e nunca compartilhados com a indústria farmacêutica ou terceiros.

        **O sistema armazena informações de pacientes?**
        Não. O serviço é uma ferramenta de atualização para o médico. Não processamos dados de seus pacientes.
        """)

    with st.expander("💳 4. GESTÃO DA ASSINATURA"):
        st.write("""
        **Posso alterar minhas keywords?**
        Sim. No Portal do Assinante, você tem autonomia total para atualizar seus focos de interesse a qualquer momento.

        **Existe fidelidade?**
        Não. Você pode cancelar sua assinatura a qualquer momento diretamente pelo painel, sem multas.
        """)

    with st.expander("📊 5. QUALIDADE DO RESUMO E DADOS"):
        st.write("""
        **A IA pode omitir dados importantes?**
        Nossa IA usa a técnica 'Chain-of-Thought', sendo obrigada a extrair Metodologia, Resultados, P-valores e Intervalos de Confiança. O objetivo é a síntese técnica fiel.

        **E se eu quiser mudar a frequência de recebimento?**
        Você está no controle. Oferecemos planos de 1x por semana a entregas diárias. Ajuste seu plano no perfil quando desejar.
        """)

    with st.expander("🌐 6. FONTES E DIRETRIZES"):
        st.write("""
        **Vocês monitoram outras bases?**
        Focamos no PubMed (padrão ouro). Também identificamos diretrizes de sociedades internacionais (AHA, ESC, AAD) indexadas, garantindo que você não perca mudanças em protocolos.
        """)

# --- ESTRUTURA DA PÁGINA ---
st.title("🩺 Portal do Assinante Medical In-Sight")
st.subheader("Gerencie sua inteligência científica personalizada")

email_login = st.sidebar.text_input("Acesse seu painel (E-mail):")

if email_login:
    dados = buscar_dados_medico(email_login)
    
    if dados:
        # Ajuste aqui: tratamos se o banco trouxer 4 ou 6 colunas
        if len(dados) == 4:
            nome, esp, keys, plano = dados
            dia, hora = "Segunda", "08:00"  # Valores padrão
        else:
            nome, esp, keys, plano, dia, hora = dados
            
        st.markdown(f"### Bem-vindo(a), Dr(a). {nome}")
        
        col1, col2, col3 = st.columns(3)
        with col1: st.metric("Seu Plano", plano)
        with col2: st.metric("Próximo Envio", f"{dia} às {hora}")
        with col3: st.metric("Especialidade", esp)

        st.divider()

        # Adicionada a TAB 3: "Histórico de Envios"
        tab1, tab2, tab3 = st.tabs(["🎯 Minha Curadoria", "🛡️ Central de Ajuda & Transparência", "📚 Histórico de Envios"])

        with tab1:
            st.write("### Ajuste seus focos de interesse")
            # Melhoria: Seleção de Especialidade agora salva no banco
            lista_especialidades = ["Cardiologia", "Psiquiatria", "Pediatria", "Neurologia", "Dermatologia", "Oncologia", "Outra"]
            try:
                indice_atual = lista_especialidades.index(esp)
            except:
                indice_atual = 0
            
            nova_esp = st.selectbox("Especialidade Principal:", lista_especialidades, index=indice_atual)
            novas_keys = st.text_input("Suas Keywords Atuais:", value=keys)
            
            if st.button("Atualizar Preferências"):
                try:
                    conn = sqlite3.connect('medical_insight.db')
                    cursor = conn.cursor()
                    cursor.execute("""
                        UPDATE clientes 
                        SET especialidade = ?, keywords = ? 
                        WHERE email = ?
                    """, (nova_esp, novas_keys, email_login))
                    conn.commit()
                    conn.close()
                    
                    st.success("✅ Preferências atualizadas!")
                    
                    # ESSE É O COMANDO QUE FAZ A PÁGINA ATUALIZAR O TOPO:
                    st.rerun() 
                    
                except Exception as e:
                    st.error(f"Erro ao salvar: {e}")

        with tab2:
            exibir_faq_completo()

        with tab3:
            st.write("### 📚 Seus Últimos Artigos Analisados")
            historico = buscar_historico(email_login)
            if historico:
                for titulo, data, link in historico:
                    # Pequeno ajuste estético para a data
                    data_f = data.split()[0] if data else ""
                    st.markdown(f"**📅 {data_f}**")
                    st.markdown(f"#### {titulo}")
                    st.link_button("📂 Abrir no PubMed", link)
                    st.divider()
            else:
                st.info("Nenhum histórico disponível ainda. Ele aparecerá após o seu primeiro envio processado.")

    else:
        st.error("E-mail não localizado. Confira os dados ou veja nossas informações abaixo.")
        st.divider()
        exibir_faq_completo()
else:
    st.info("Digite seu e-mail na barra lateral para acessar seu histórico. Se você ainda não é assinante, conheça nosso compromisso abaixo:")
    st.divider()
    exibir_faq_completo()

# Rodapé
st.markdown("---")
st.markdown("<p style='text-align: center; color: gray;'>Medical In-Sight © 2026 | Protocolo v2.4.0 | Tecnologia de IA validada por profissionais reais.</p>", unsafe_allow_html=True)