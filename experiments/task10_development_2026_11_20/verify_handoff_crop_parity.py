"""Existing development VCF crop and preparation parity; no genotype resimulation."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
from structural_inference_handoff import crop_to_map
from scripts.develop_task10_detected_ibd import save_json, sha256, verify_inputs
from scripts.prepare_hapibd_vcf import prepare_vcf


def content_hash(path):
    digest=hashlib.sha256()
    with gzip.open(path,'rb') as stream:
        while data:=stream.read(1<<20):
            digest.update(data)
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Scheduled parity job required')
    protocol=json.loads((args.source/'structural_inference_handoff_protocol.json').read_text())
    manifest=json.loads((args.source/'source_manifest.json').read_text())
    verify_inputs(manifest)
    output=args.root/'handoff_crop_parity_v1'
    output.mkdir(parents=True,exist_ok=False)
    reports=[]
    bound={}
    for index in range(40,48):
        folder=args.root/'detected_panel_v2'/f'cohort_{index:03d}'
        original=folder/'family.normalized.vcf.gz'
        prepared=folder/'prepared.vcf.gz'
        bound[str(original)]=sha256(original)
        bound[str(prepared)]=sha256(prepared)
        cropped=output/f'crop_{index}.vcf.gz'
        destination=output/f'prepared_{index}.vcf.gz'
        crop=crop_to_map(original,cropped,protocol['map_bounds'],protocol['minimum_markers_per_chromosome'])
        contract=prepare_vcf(cropped,destination,missing_policy='error')
        comparison={'original_content_sha256':content_hash(original),'cropped_content_sha256':content_hash(cropped),'old_prepared_content_sha256':content_hash(prepared),'new_prepared_content_sha256':content_hash(destination)}
        if comparison['original_content_sha256']!=comparison['cropped_content_sha256'] or comparison['old_prepared_content_sha256']!=comparison['new_prepared_content_sha256'] or crop['records_excluded']:
            raise ValueError('Crop or preparation altered existing qualification-domain inputs')
        reports.append({'index':index,'crop':crop,'contract':contract,'comparison':comparison})
    verify_inputs(bound)
    save_json(output/'result.json',{'state':'PASS','cohorts':reports,'input_sha256':bound,'source_manifest_sha256':sha256(args.source/'source_manifest.json'),'protocol_sha256':sha256(args.source/'structural_inference_handoff_protocol.json'),'consequence':'Existing detector inputs byte-identical after decompression, therefore no changed inference input; no new biological evidence'})


if __name__=='__main__':
    main()
