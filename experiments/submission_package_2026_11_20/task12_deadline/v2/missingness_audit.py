from pathlib import Path
import gzip,json,hashlib,datetime
p=Path('/private/tmp/ghist-final-submissions-20261120/task12/relatedness4.final.vcf.gz')
x={'records':0,'missing_sites':0,'missing_alleles':0,'genotype_alleles':0,'unphased_heterozygotes':0,'biallelic_sites':0,'biallelic_missing_sites':0,'all_missing_sites':0,'missing_by_sample':{}}
with gzip.open(p,'rt') as h:
 for line in h:
  if line.startswith('#CHROM'):
   samples=line.rstrip().split('\t')[9:];x['samples']=len(samples)
  if line.startswith('#'):continue
  f=line.rstrip().split('\t');x['records']+=1
  biallelic=',' not in f[4];x['biallelic_sites']+=biallelic
  idx=f[8].split(':').index('GT');missing=0;called=0
  for name,cell in zip(samples,f[9:],strict=True):
   gt=cell.split(':')[idx];alleles=gt.replace('|','/').split('/');n=alleles.count('.');missing+=n;called+=len(alleles)-n;x['genotype_alleles']+=len(alleles)
   if n:x['missing_by_sample'][name]=x['missing_by_sample'].get(name,0)+n
   if '/' in gt and len(alleles)==2 and '.' not in alleles and alleles[0]!=alleles[1]:x['unphased_heterozygotes']+=1
  x['missing_alleles']+=missing;x['missing_sites']+=missing>0;x['biallelic_missing_sites']+=biallelic and missing>0;x['all_missing_sites']+=called==0
x['input_sha256']=hashlib.sha256(p.read_bytes()).hexdigest();x['created_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();x['scope']='Input genotype contract measurement only; no phenotype/truth/relationship inference.'
out=p.with_name('missingness_audit.json');out.write_text(json.dumps(x,indent=2)+'\n');print(json.dumps({k:v for k,v in x.items() if k!='missing_by_sample'},indent=2))
