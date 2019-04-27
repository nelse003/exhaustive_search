import sqlite3
import pandas as pd

# NUDT22
# database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"

# NUDT7
# database_path = "/dls/labxchem/data/2017/lb18145-49/processing/database/soakDBDataFile.sqlite"

# NUDT7 (Tobias)
# database_path = "/dls/labxchem/data/2017/lb18145-3/processing/database/soakDBDataFile.sqlite"

# FALZA
database_path = (
    "/dls/labxchem/data/2016/lb13385-61/processing/database/soakDBDataFile.sqlite"
)

conn = sqlite3.connect(database_path)
main_table_df = pd.read_sql_query("select * from mainTable", conn)
cur = conn.cursor()

refinement_outcomes = (
    " '3 - In Refinement', '4 - CompChem ready', '5 - Deposition ready','6 - Deposited'"
)

cur.execute(
    "SELECT CrystalName, CompoundCode, RefinementResolution "
    "FROM mainTable WHERE RefinementOutcome in ({})"
    " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(
        refinement_outcomes
    )
)

refinement_xtals = cur.fetchall()
cur.close()

compounds = {}
for xtal_name, compound_code, resolution in refinement_xtals:
    compounds[xtal_name] = compound_code

comp_df = pd.DataFrame(
    list(compounds.items()), columns=["CrystalName", "compound_code"]
)
duplicate_compound_df = pd.concat(
    g for _, g in comp_df.groupby("compound_code") if len(g) > 1
)

print(duplicate_compound_df["compound_code"].value_counts())
