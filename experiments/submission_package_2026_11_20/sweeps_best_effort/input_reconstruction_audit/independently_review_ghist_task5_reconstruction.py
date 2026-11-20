"""Read-only paired audit of Task5 coordinate reconstruction; no inference."""
from pathlib import Path
import datetime
import gzip
import hashlib
import itertools
import json
import re
import pysam

base = Path('/private/tmp/ghist-final-submissions-20261120/task5')
source = base/'onesweep.final.original_v1.vcf.gz'
target = base/'onesweep.final.reconstructed_pos1.vcf.gz'
producer = Path('/tmp/reconstruct_ghist_task5_pos1.py')
provenance_path = base/'reconstructed_pos1_provenance.json'

def digest(path, algorithm='sha256'):
    h=hashlib.new(algorithm)
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1<<20),b''):h.update(block)
    return h.hexdigest()

pre={str(p):digest(p) for p in (source,target,producer,provenance_path)}
provenance=json.loads(provenance_path.read_text())
assert pre[str(source)]=='4ffd4b17fde78c182fd79b73ea7ebb246acb077d0bdc75fc04ff65c905377a14'
assert digest(source,'md5')=='242445cad013f18cf7c661a869ccd2d1'
assert pre[str(target)]==provenance['target_sha256']
assert pre[str(producer)]==provenance['script_sha256']
assert provenance['current_release_bytes_verified'] is False
assert provenance['current_release']['release_note']=='Correct VCF POS to 1-based (+1); all other columns unchanged'
headers=hashlib.sha256(); old_content=hashlib.sha256();new_content=hashlib.sha256()
record_count=header_count=0
last_old=last_new=None;first_old=first_new=None
lengths={};sample_count=None
with gzip.open(source,'rb') as original,gzip.open(target,'rb') as derived:
    for line_number,(old,new) in enumerate(itertools.zip_longest(original,derived),1):
        assert old is not None and new is not None,('unequal line count',line_number)
        old_content.update(old);new_content.update(new)
        if old.startswith(b'#'):
            assert old==new,('header changed',line_number)
            assert record_count==0,'header after records'
            header_count+=1;headers.update(old)
            if old.startswith(b'##contig='):
                match=re.fullmatch(rb'##contig=<ID=([^,>]+),length=(\d+)>\r?\n',old)
                assert match,('unexpected contig syntax',old)
                assert match[1] not in lengths,'duplicate contig'
                lengths[match[1]]=int(match[2])
            if old.startswith(b'#CHROM'):
                sample_count=len(old.rstrip(b'\r\n').split(b'\t'))-9
            continue
        assert not new.startswith(b'#')
        left=old.split(b'\t');right=new.split(b'\t')
        assert len(left)==len(right)==9+sample_count,('column count',line_number)
        assert left[0]==right[0]==b'9'
        assert left[2:]==right[2:],('non-POS field changed',line_number)
        a,b=int(left[1]),int(right[1]);assert b==a+1,('POS transformation',line_number)
        assert 0<=a<lengths[b'9'] and 1<=b<=lengths[b'9']
        assert last_old is None or a>last_old,'original duplicate or unsorted position'
        assert last_new is None or b>last_new,'derived duplicate or unsorted position'
        if first_old is None:first_old,first_new=a,b
        last_old,last_new=a,b;record_count+=1
assert lengths=={b'9':138394717}
assert sample_count==25 and record_count>0
assert record_count==provenance['records']==provenance['paired_records_verified']
assert [first_old,last_old]==provenance['original_position_range']
assert [first_new,last_new]==provenance['corrected_position_range']
try:
    with pysam.VariantFile(str(target)) as vcf:
        next(iter(vcf))
    parser={'compatible':True,'pysam_version':pysam.__version__}
except Exception as error:
    parser={'compatible':False,'pysam_version':pysam.__version__,'error_type':type(error).__name__,'error':str(error)}
assert {str(p):digest(p) for p in (source,target,producer,provenance_path)}==pre,'file changed during audit'
report={
 'schema':'independent_task5_pos_reconstruction_review_v1',
 'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'reviewer':'sweeps_science independent paired read-only audit',
 'review_script_path':str(Path(__file__)), 'review_script_sha256':digest(Path(__file__)),
 'bound_files_sha256':pre,
 'source_entity':'syn75346852','source_version':1,
 'source_live_metadata_attestation':'Parent independently authenticated original download and supplied live MD5; this audit rechecks local MD5/SHA256, not browser state.',
 'source_md5':digest(source,'md5'), 'derived_md5':digest(target,'md5'),
 'header_lines':header_count,'headers_identical':True,'header_bytes_sha256':headers.hexdigest(),
 'records':record_count,'paired_record_count_equal':True,'all_non_POS_fields_byte_identical':True,
 'all_POS_exactly_plus_one':True,'strictly_increasing_unique_positions':True,
 'original_position_range':[first_old,last_old],'derived_position_range':[first_new,last_new],
 'contig_lengths':{'9':138394717},'diploid_sample_columns':sample_count,
 'source_uncompressed_sha256':old_content.hexdigest(),'derived_uncompressed_sha256':new_content.hexdigest(),
 'scientific_parameters_changed':False,'input_content_review':'PASS: exact documented POS-only transformation for every row',
 'locked_pysam_reader_smoke':parser,
 'operational_verdict':'BLOCKED for direct use of ordinary gzip artifact; retain it, derive plain VCF or BGZF with identical uncompressed bytes, record separate hashes, and verify parser compatibility before inference.' if not parser['compatible'] else 'Parser compatible',
 'current_release_version':2,'current_release_bytes_verified':False,
 'current_release_content_equivalence':'INFERRED ONLY from parent-supplied official release note: Correct VCF POS to 1-based (+1); all other columns unchanged. This audit has no current-v2 bytes and does not prove release identity.',
 'public_answer_or_truth_access':False,'inference_or_jobs_launched':False,'external_submission_performed':False,
}
output=base/'independent_reconstruction_review.json'
assert not output.exists()
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'review_path':str(output),'review_sha256':digest(output),'records':record_count,'headers':header_count,'position_range':[first_new,last_new],'parser':parser,'current_release_bytes_verified':False},indent=2))
