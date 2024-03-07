#!/usr/bin/env -S sage -python

"""
This script is used to convert perform additional verification on uploads before review by a human editor.  Simple verifications are performed when an LMFDB user uploads data, but additional verification (anything requiring a nontrivial amount of computation for example) are performed by this script.  It is intended for use by a cron job, and is called with no arguments.
"""

import os
import sys
import tempfile
from collections import defaultdict
from datetime import datetime
here = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(here, "data")
upone, _ = os.path.split(here)
uptwo, _ = os.path.split(upone)
sys.path.append(uptwo)
from lmfdb import db
from lmfdb.backend.encoding import copy_dumps
from lmfdb.users.main import Reviewer


def process_all():
    os.makedirs("data", exist_ok=True)
    by_table = defaultdict(list)
    cols = {}
    ids = []
    reviewer = Reviewer()
    # TODO: it might be better to isolate each verification in its own process to insulate against timeouts
    # This would prevent sharing some computations, like the poset of modular curves, between different uploads
    with tempfile.NamedTemporaryFile("w", delete=False) as F:
        columns = ["id", "status", "processed", "updated", "comment"]
        types = ["bigint", "smallint", "timestamp without time zone", "timestamp without time zone", "text"]
        _ = F.write("|".join(columns) + "\n" + "|".join(types) + "\n\n")
        for rec in db.data_uploads.search({"status":2}, ["data", "id", "section"]):
            try:
                section = reviewer.section_lookup[rec["section"]]
                for table, newrow, line in section.process(rec["data"]):
                    if (table, newrow) in cols:
                        if cols[table, newrow] != set(line):
                            raise ValueError(f"Schema for id {rec['id']}:{table}:{newrow} does not agree with previous schema")
                    else:
                        cols[table, newrow] = set(line)
                    col_type = db[table].col_type
                    line = "|".join(copy_dumps(line[col], col_type[col]) for col in sorted(line))
                    by_table[table, newrow].append(line)
            except Exception as err:
                status = -3
                comment = str(err)
            else:
                status = 3
                comment = ""
                ids.append(rec["id"])
            timestamp = datetime.utcnow().isoformat()
            _ = F.write(f"{rec['id']}|{status}|{timestamp}|{timestamp}|{comment}\n")
        F.close()
        db.data_uploads.update_from_file(F.name, "id")
        db.data_uploads.cleanup_from_reload()
        os.unlink(F.name)
    timestamp = datetime.utcnow().isoformat().replace(":", "-").replace("T", "-").replace(".", "-")
    uploads = []
    for (table, newrows), lines in by_table.items():
        nr = "t" if newrows else "f"
        fname = os.path.join("data", f"{table}_{nr}_{timestamp}")
        columns = sorted(cols[table, newrows])
        col_type = db[table].col_type
        types = [col_type[col] for col in columns]
        with open(fname, "w") as F:
            _ = F.write("|".join(columns) + "\n" + "|".join(types) + "\n\n")
            for line in lines:
                _ = F.write(line + "\n")
        uploads.append((table, newrows, fname))
    try:
        for table, newrows, fname in uploads:
            if newrows:
                db[table].copy_from(fname)
            else:
                db[table].update_from_file(fname)
    except Exception as err:
        payload = {"status": -4, "comment": str(err)}
    else:
        payload = {"status": 4}
    db.data_uploads.update({"id": {"$in": ids}}, payload)

process_all()
