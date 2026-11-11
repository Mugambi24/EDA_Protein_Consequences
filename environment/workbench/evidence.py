"""Independent RNA junction-support audit module (buggy baseline)."""

def summarize(data):
    """Aggregate per-sample molecules with quality filtering, UMI deduplication,
    prior-smoothed A1 proportions, uncertainty and case/control comparisons.

    This incomplete baseline intentionally counts raw reads and must be repaired.
    """
    table={}
    for row in data['observations']:
        cid=row['case_id']
        if cid not in table: table[cid]={'case_id':cid,'raw_count':0}
        table[cid]['raw_count']+=1
    return [table[cid] for cid in sorted(table)]
