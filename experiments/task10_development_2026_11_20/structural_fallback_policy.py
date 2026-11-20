"""Explicit pragmatic whole-cohort fallback for one declared inference condition."""
from __future__ import annotations
import hashlib
import itertools
import json
from structural_pair_candidate import infer as structural_infer

NO_ROOT = 'Initial NJ tree has no accepted root; candidate cannot provide structural predictions'


def canonical_complete(rows):
    samples=sorted({str(row[key]) for row in rows for key in ['indID1','indID2']})
    expected=set(itertools.combinations(samples,2))
    normalized=[]
    for row in rows:
        a,b=sorted([row['indID1'],row['indID2']])
        normalized.append(dict(row,indID1=a,indID2=b))
    normalized.sort(key=lambda row:(row['indID1'],row['indID2']))
    observed=[(row['indID1'],row['indID2']) for row in normalized]
    if len(samples)<2 or len(observed)!=len(expected) or set(observed)!=expected:
        raise ValueError('Pragmatic policy requires complete unique unordered roster pairs')
    return normalized,samples


def infer_with_fallback(pair_rows):
    canonical,samples=canonical_complete(pair_rows)
    provenance={'policy':'structural-if-valid-else-exact-initial-no-root-whole-cohort-pair-only-v1',
        'people':len(samples),'pairs':len(canonical),'input_pair_rows_sha256':hashlib.sha256(json.dumps(canonical,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        'sample_order':samples,'strict_qualification_claim':False,'individual_uncertainty_calibrated':False}
    try:
        rows,tree=structural_infer(canonical)
    except ValueError as error:
        if str(error)!=NO_ROOT:
            raise
        provenance.update({'used_fallback':True,'selected_method':'frozen_pair_only','reason':str(error),
            'selection_unit':'entire cohort; no truth or pairwise selective replacement'})
        return canonical,{'pragmatic_policy':provenance,'structural_tree_unavailable':True}
    ordered,observed=canonical_complete(rows)
    if observed!=samples:
        raise ValueError('Structural output roster changed')
    provenance.update({'used_fallback':False,'selected_method':'frozen_structural','reason':None})
    return ordered,dict(tree,pragmatic_policy=provenance)
