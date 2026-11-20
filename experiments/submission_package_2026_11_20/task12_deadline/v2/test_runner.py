import importlib.util
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('task12final',Path(__file__).with_name('run.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class Contracts(unittest.TestCase):
 def test_zero_pairs_preserve_present_correction(self):
  from scripts.fit_human_family_relationships import summarize_raw_lengths
  from scripts.fit_human_background_corrected import correct_records
  lock={'genome_cm':3346.298892538,'minimum_cm':5}
  ref=[{'indID1':'r1','indID2':'r2',**summarize_raw_lengths([('1',8)],lock['genome_cm'])}]
  grouped={('a','b'):[('1',8),('1',12)]}
  rows,n=m.feature_records(grouped,['c','a','b'],ref,lock)
  original=[{'indID1':'a','indID2':'b',**summarize_raw_lengths(grouped[('a','b')],lock['genome_cm'])}]
  self.assertEqual(rows[0],correct_records(original,ref,5,lock['genome_cm'])[0]);self.assertEqual(n,2)
  self.assertTrue(all(v==0 for r in rows[1:] for k,v in r.items() if k not in {'indID1','indID2'}))
 def test_all_zero_pairs_do_not_quantile_empty(self):
  rows,n=m.feature_records({},['a','b','c'],[],{'genome_cm':3346.298892538,'minimum_cm':5})
  self.assertEqual((len(rows),n),(3,3))
 def test_probability_and_coverage_guard(self):
  from scripts.fit_human_family_relationships import LABELS
  p=dict.fromkeys(LABELS,0.);p[LABELS[-1]]=1.
  row={'indID1':'a','indID2':'b','relation':LABELS[-1],'confidence':1.,'probabilities':p}
  m.check_predictions([row],['a','b'])
  with self.assertRaises(ValueError):m.check_predictions([row],['a','b','c'])
  p[LABELS[0]]=float('nan')
  with self.assertRaises(ValueError):m.check_predictions([row],['a','b'])
 def test_explicit_major_preserves_phase_and_counts(self):
  from scripts.prepare_hapibd_vcf import prepare_vcf
  with tempfile.TemporaryDirectory() as tmp:
   source=Path(tmp)/'x.vcf';out=Path(tmp)/'y.vcf'
   source.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t.|.\n')
   with self.assertRaises(ValueError):prepare_vcf(source,out,missing_policy='error')
   self.assertFalse(out.exists())
   report=prepare_vcf(source,out,missing_policy='major')
   self.assertEqual(report['imputed_alleles'],2)
   self.assertTrue(out.read_text().endswith('0|1\t0|0\n'))
 def test_protocol_binding_guard(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'protocol.json';p.write_text('{"seed":1}')
   _,h=m.read_bound_json(p);p.write_text('{"seed":2}')
   self.assertNotEqual(h,m.digest(p))
 def test_vcf_coordinates(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'x.vcf'
   header='##fileformat=VCFv4.2\n##contig=<ID=1,length=20>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n'
   p.write_text(header+'1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n')
   self.assertEqual(m.inspect_vcf(p)[0],['a','b'])
   p.write_text(header+'1\t0\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n')
   with self.assertRaises(ValueError):m.inspect_vcf(p)

if __name__=='__main__':unittest.main()
