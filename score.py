#!/usr/bin/env python3
"""Score search results against the SCOP labels carried by the released data.

Labels come straight from the data directory:
  * database.fasta: annotated sequences look like
        >annotated_C9ST65 g.78.1.1 |1.7e-06 |g1sse.1 |35-121 |349-433
    (a full-length protein may carry several SCOP domains before the first
    "|"); every other sequence is a reversed decoy and has no label.
  * query/*.fasta: the header carries the query's own SCOP code(s).

Judging convention:

  hit                                          verdict
  -------------------------------------------- ------------------------------
  any domain in the query's superfamily         TP  (counts toward recall)
  otherwise any domain in the query's fold      ignored (homology disputed)
  otherwise query or hit has a class e domain   ignored (see below)
  otherwise labeled (different fold)            FP
  no label at all (reversed decoy)              FP

SCOP class e ("multi-domain proteins") groups domains whose assignment is
ambiguous, so a non-TP match involving it is not counted as a false positive.
This rule is taken from the MMseqs2 benchmark evaluator (EvaluateResults.cpp).

Metrics per query, then averaged over all official queries:
  recall    = TP / (database sequences in the query's superfamily)
  precision = TP / (TP + FP)
  ROC1      = true homologs ranked above the first FP, / all true homologs
  ROC10     = Gribskov-Robinson ROC_n, n=10: mean over the first 10 FPs of
              "true homologs ranked above it", / all true homologs.
              Fewer than 10 FPs: the curve stays flat at the final TP count.
ROC metrics READ THE HIT FILE IN ORDER -- write hits most significant first.
ROC10 is the primary quality metric.

Every file in query/ is an official query. A query without a hit table
scores 0 on every metric.

Usage: score.py [--data DIR] [--results DIR] [--format tbl|m8] [-o FILE]
  --data     released data directory (default: ../preliminary next to the repo)
  --results  run_queries.sh output directory (default: results/ in the repo)
Reads <results>/tbl/*.tbl (hmmer --tblout) or *.m8 (query,target in cols 1,2).
Official scoring is recomputed offline by the organizers; submitted hit tables
must be produced by your actual search program.
"""
import argparse
import collections
import glob
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.abspath(__file__))
SCCS = re.compile(r'^[a-z]\.\d+\.\d+\.\d+$')
ROC_N = 10


def superfamily(code):
    return code.rsplit('.', 1)[0]


def fold(code):
    return '.'.join(code.split('.')[:2])


def header_codes(fields):
    """SCOP codes in a header, i.e. the tokens before the first '|' field."""
    codes = []
    for tok in fields:
        if tok.startswith('|'):
            break
        if SCCS.match(tok):
            codes.append(tok)
    return codes


def annotated_headers(db):
    """Yield the header lines of labeled sequences (grep when available)."""
    try:
        p = subprocess.Popen(['grep', '^>annotated_', db], stdout=subprocess.PIPE,
                             env=dict(os.environ, LC_ALL='C'))
    except OSError:
        with open(db) as fh:
            for line in fh:
                if line.startswith('>annotated_'):
                    yield line
        return
    for line in p.stdout:
        yield line.decode('ascii', 'replace')
    p.wait()


def load_labels(db):
    """-> id -> (superfamilies, folds, has class e); superfamily -> member count."""
    labels = {}
    sf_size = collections.Counter()
    for line in annotated_headers(db):
        f = line.split()
        codes = header_codes(f[1:])
        if not codes:
            continue
        sfs = frozenset(superfamily(c) for c in codes)
        labels[f[0][1:]] = (sfs, frozenset(fold(c) for c in codes),
                            any(c.startswith('e.') for c in codes))
        for s in sfs:
            sf_size[s] += 1
    return labels, sf_size


def load_queries(qdir):
    """-> query id -> (superfamilies, folds, has class e), from the query headers."""
    queries = {}
    for path in sorted(glob.glob(os.path.join(qdir, '*.fasta'))):
        with open(path) as fh:
            f = fh.readline()[1:].split()
        codes = header_codes(f[1:])
        if codes:
            queries[f[0]] = (frozenset(superfamily(c) for c in codes),
                             frozenset(fold(c) for c in codes),
                             any(c.startswith('e.') for c in codes))
    return queries


def read_hits(paths, fmt):
    """-> query -> target ids in file (rank) order, deduplicated."""
    hits = collections.defaultdict(list)
    seen = collections.defaultdict(set)
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                f = line.split()
                q, t = (f[0], f[1]) if fmt == 'm8' else (f[2], f[0])
                if t not in seen[q]:
                    seen[q].add(t)
                    hits[q].append(t)
    return hits


def score_query(q, qsf, qfold, q_e, targets, labels, sf_size):
    n_true = sum(sf_size[s] for s in qsf)
    tp = fp = ign = 0
    first_fp = 0
    tp_above_fp = []          # true homologs ranked above the i-th FP
    for rank, t in enumerate(targets, 1):
        if t == q:
            continue
        lab = labels.get(t)
        if lab and lab[0] & qsf:
            tp += 1
            continue
        if lab and (lab[1] & qfold or lab[2] or q_e):
            ign += 1
            continue
        fp += 1
        if not first_fp:
            first_fp = rank
        if len(tp_above_fp) < ROC_N:
            tp_above_fp.append(tp)
    while len(tp_above_fp) < ROC_N:
        tp_above_fp.append(tp)
    recall = tp / n_true if n_true else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    roc1 = tp_above_fp[0] / n_true if n_true else 0.0
    roc10 = sum(tp_above_fp) / (ROC_N * n_true) if n_true else 0.0
    return (q, ','.join(sorted(qsf)), n_true, len(targets), tp, fp, ign,
            first_fp, recall, precision, roc1, roc10)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.join(os.path.dirname(REPO), 'preliminary'),
                    help='released data directory (query/ + database.fasta)')
    ap.add_argument('--results', default=os.path.join(REPO, 'results'),
                    help='run_queries.sh output directory')
    ap.add_argument('--format', choices=['tbl', 'm8'], default='tbl')
    ap.add_argument('-o', '--output', default=None)
    args = ap.parse_args()

    db = os.path.join(args.data, 'database.fasta')
    qdir = os.path.join(args.data, 'query')
    for p in (db, qdir):
        if not os.path.exists(p):
            sys.exit('ERROR: %s not found (set --data)' % p)
    ext = '*.m8' if args.format == 'm8' else '*.tbl'
    files = sorted(glob.glob(os.path.join(args.results, 'tbl', ext)))
    if not files:
        sys.exit('ERROR: no %s files under %s' % (ext, os.path.join(args.results, 'tbl')))

    queries = load_queries(qdir)
    print('reading labels from %s ...' % db, flush=True)
    labels, sf_size = load_labels(db)
    hits = read_hits(files, args.format)

    # A query "has a hit table" if one was written for it, even with zero hits
    # (run_queries.sh names each file <query>.tbl).
    produced = set(hits) | set(os.path.basename(p).rsplit('.', 1)[0] for p in files)
    rows = []
    missing = []
    for q in sorted(queries):
        if q not in produced:
            missing.append(q)
        qsf, qfold, q_e = queries[q]
        rows.append(score_query(q, qsf, qfold, q_e, hits.get(q, []), labels, sf_size))
    extra = sorted(set(hits) - set(queries))

    out_path = args.output or os.path.join(args.results, 'scores.tsv')
    with open(out_path, 'w') as fh:
        fh.write('query\tsuperfamily\ttrue_homologs\thits\ttp\tfp\tignored'
                 '\tfirst_fp_rank\trecall\tprecision\troc1\troc10\n')
        for r in rows:
            fh.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n' % r)

    n = len(rows)
    mean = lambda i: sum(r[i] for r in rows) / n
    print('labeled sequences in database: %d  (%d superfamilies)' % (len(labels), len(sf_size)))
    print('official queries: %d  |  with hit table: %d  |  detail: %s'
          % (n, n - len(missing), out_path))
    if missing:
        print('WARNING: %d official queries have no hit table and score 0: %s%s'
              % (len(missing), ' '.join(missing[:5]), ' ...' if len(missing) > 5 else ''))
    if extra:
        print('note: %d hit-table queries are not official and were ignored' % len(extra))
    print('mean recall    = %.4f' % mean(8))
    print('mean precision = %.4f' % mean(9))
    print('mean ROC1      = %.4f' % mean(10))
    print('mean ROC10     = %.4f   <- primary quality metric' % mean(11))

    wall = os.path.join(args.results, 'wall_time.tsv')
    if os.path.exists(wall):
        with open(wall) as fh:
            fh.readline()
            f = fh.readline().rstrip('\n').split('\t')
        try:
            secs, ran, skipped = float(f[0]), int(f[1]), int(f[2])
        except (IndexError, ValueError):
            print('WARNING: cannot parse %s' % wall)
        else:
            print('search time    = %.1f s   <- official: wall clock from end of warm-up '
                  'to last query (%d run)' % (secs, ran))
            if skipped:
                print('WARNING: %d queries reused earlier results; this time is NOT a valid '
                      'official time (clear the output dir or set FORCE=1)' % skipped)
    else:
        print('search time    = n/a (no wall_time.tsv in %s)' % args.results)

if __name__ == '__main__':
    main()
