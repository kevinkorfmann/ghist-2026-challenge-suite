import collections,json,os,pathlib,subprocess
import numpy as np,pysam,tskit
root=pathlib.Path('/vast/projects/smathi/cohort/korfmann/ghist2026')
out=root/'scientific/readiness_sweeps/export_contract';out.mkdir(exist_ok=True)
ts=tskit.load(root/'artifacts/task5_single_sweep/selection_exact_recovery/task_0.trees')
p=out/'original.vcf'
with p.open('w') as f: ts.write_vcf(f,contig_id='7',position_transform=lambda x:np.floor(x).astype(int)+1)
labels=collections.Counter(); biallelic=usable=records=0
with pysam.VariantFile(p) as f:
 length=f.header.contigs['7'].length; samples=len(f.header.samples)
 for row in f:
  records+=1; labels.update(row.alleles)
  if len(row.alleles)==2:
   biallelic+=1
   calls=[v for s in row.samples.values() for v in s['GT']]
   if all(v in (0,1) for v in calls) and 0<sum(calls)<len(calls): usable+=1
with (out/'raisd.log').open('w') as log:
 result=subprocess.run([str(root/'external/raisd/RAiSD'),'-n','contract','-I',str(p),'-w','50','-R','-f'],cwd=out,stdout=log,stderr=subprocess.STDOUT)
output={'tskit':tskit.__version__,'sequence_length':ts.sequence_length,'vcf_contig_length':length,'samples':samples,'records':records,'biallelic_records':biallelic,'usable_biallelic_polymorphic_records':usable,'allele_label_counts':dict(labels),'raisd_returncode':result.returncode,'raisd_log':(out/'raisd.log').read_text()}
(out/'contract.json').write_text(json.dumps(output,indent=2)+'\n'); print(json.dumps(output,indent=2))
