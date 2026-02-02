import os
import sqlite3
import psycopg2
from dotenv import load_dotenv

load_dotenv()

def get_connection():
    """
    Gerencia a conexão de forma híbrida: 
    1. Se houver DATABASE_URL no .env, conecta ao PostgreSQL (Nuvem).
    2. Caso contrário, usa o arquivo medical_insight.db (Local).
    """
    database_url = os.getenv("DATABASE_URL")
    
    if database_url:
        # Conexão para Nuvem (PostgreSQL)
        # O sslmode='require' é obrigatório na maioria dos serviços cloud (Supabase, Railway)
        return psycopg2.connect(database_url, sslmode='require')
    else:
        # Conexão para Desenvolvimento Local (SQLite)
        # Mantém a compatibilidade com o que você já construiu
        return sqlite3.connect('medical_insight.db')