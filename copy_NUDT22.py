import os
import sqlite3

out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT22_repeat_soaks"
database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"
refinement_outcomes= " '3 - In Refinement', '4 - CompChem ready', '5 - Deposition ready','6 - Deposited'"

conn = sqlite3.connect(database_path)
cur = conn.cursor()

cur.execute("SELECT CrystalName, RefinementPDB_latest, "
            "RefinementMTZ_latest, CompoundCode "
            "FROM mainTable WHERE (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL")

refinement_xtals = cur.fetchall()

for xtal_name, pdb, mtz, compound_code in refinement_xtals:

    compound_dir = os.path.join(out_path, compound_code)

    if not os.path.exists(compound_dir):
        os.mkdir(compound_dir)

    xtal_dir = os.path.join(compound_dir, xtal_name)

    if not os.path.exists(xtal_dir):
        os.mkdir(xtal_dir)

    os.symlink(pdb, os.path.join(xtal_dir, "refine.pdb"))

    os.symlink(mtz, os.path.join(xtal_dir, "refine.mtz"))


