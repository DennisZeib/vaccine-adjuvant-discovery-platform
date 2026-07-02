import sqlite3
import pandas as pd
from datetime import datetime

DB_FILE = "simulations.db"

def init_db():
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS simulations (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp TEXT,
            microbe_a TEXT,
            microbe_b TEXT,
            glucose INTEGER,
            duration INTEGER,
            ratio_a INTEGER,
            ratio_b INTEGER,
            lps REAL,
            mannan REAL,
            flagellin REAL,
            ethanol REAL,
            lactate REAL,
            acetate REAL,
            score REAL,
            toxicity REAL
        )
    ''')
    conn.commit()
    conn.close()

def save_simulation(microbe_a, microbe_b, glucose, duration, ratio_a, ratio_b, final, score, toxicity):
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute('''
        INSERT INTO simulations (
            timestamp, microbe_a, microbe_b, glucose, duration,
            ratio_a, ratio_b, lps, mannan, flagellin,
            ethanol, lactate, acetate, score, toxicity
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (
        datetime.now().isoformat(),
        microbe_a, microbe_b, glucose, duration,
        ratio_a, ratio_b,
        final.get("polysaccharide_lps", 0),
        final.get("mannan_polysaccharide", 0),
        final.get("flagellin", 0),
        final.get("ethanol", 0),
        final.get("lactate", 0),
        final.get("acetate", 0),
        score,
        toxicity
    ))
    conn.commit()
    conn.close()

def load_history():
    conn = sqlite3.connect(DB_FILE)
    df = pd.read_sql_query("SELECT * FROM simulations ORDER BY id DESC", conn)
    conn.close()
    return df

def delete_all():
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute("DELETE FROM simulations")
    conn.commit()
    conn.close()