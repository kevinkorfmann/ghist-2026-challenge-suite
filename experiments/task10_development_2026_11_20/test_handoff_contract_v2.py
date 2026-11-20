import gzip
import itertools
import pytest
from structural_inference_handoff_v2 import crop_to_map, write_candidate_csv
from run_site_thinning_diagnostic import thin


def vcf(path, positions, chromosome='1'):
    with gzip.open(path,'wt') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n')
        for pos in positions:
            f.write(f'{chromosome}\t{pos}\t.\tA\tG\t.\t.\t.\tGT\t0|1\t1|0\n')


def test_crop_includes_endpoints_preserves_alleles_and_identity(tmp_path):
    a,b,c= [tmp_path / name for name in ['a.gz','b.gz','c.gz']]
    vcf(a,[9,10,15,20,21])
    report=crop_to_map(a,b,{'1':[10,20]},3)
    assert report['records_excluded']==2
    assert [int(x.split('\t')[1]) for x in gzip.open(b,'rt') if not x.startswith('#')]==[10,15,20]
    crop_to_map(b,c,{'1':[10,20]},3)
    assert gzip.open(b,'rb').read()==gzip.open(c,'rb').read()
    assert '0|1\t1|0' in gzip.open(c,'rt').read()


def test_invalid_crop_never_publishes(tmp_path):
    a,b=tmp_path/'a.gz',tmp_path/'b.gz'
    vcf(a,[10,9])
    with pytest.raises(ValueError,match='increasing'):
        crop_to_map(a,b,{'1':[1,20]},1)
    assert not b.exists()
    vcf(a,[10])
    with pytest.raises(ValueError,match='inadequate'):
        crop_to_map(a,b,{'1':[1,20],'2':[1,20]},1)
    assert not b.exists()


@pytest.mark.parametrize("people", [229, 263])
def test_csv_all_original_ids_and_missing_pair_rejection(tmp_path, people):
    names=[f'original-{i:03}' for i in range(people)]
    rows=[{'indID1':a,'indID2':b,'relation':'siblings'} for a,b in itertools.combinations(names,2)]
    out=tmp_path/'candidate.csv'
    write_candidate_csv(rows,names,out)
    assert len(out.read_text().splitlines())==1+people*(people-1)//2
    with pytest.raises(ValueError,match='every unique'):
        write_candidate_csv(rows[:-1],names,tmp_path/'bad.csv')
    assert not (tmp_path/'bad.csv').exists()


def test_thinning_deterministic_exact_rate_and_no_call_changes(tmp_path):
    a,b,c= [tmp_path / name for name in ['a.gz','b.gz','c.gz']]
    vcf(a,range(1,101))
    report=thin(a,b,40)
    assert report['retained']=={'1':75}
    assert thin(a,c,40)==report
    assert gzip.open(b,'rb').read()==gzip.open(c,'rb').read()
    assert set(gzip.open(b,'rt')).issubset(set(gzip.open(a,'rt')))
