import os
import pandas as pd
import sqlite3

database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"

conn = sqlite3.connect(database_path)
main_table_df = pd.read_sql_query("select * from mainTable",conn)
cur = conn.cursor()

cur.execute("SELECT CrystalName, CompoundCode, RefinementResolution "
            "FROM mainTable WHERE RefinementOutcome in ({})"
            " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(refinement_outcomes))

refinement_xtals = cur.fetchall()
cur.close()

compounds = {}
for xtal_name, compound_code, resolution in refinement_xtals:
    compounds[xtal_name] = compound_code

print(compounds)