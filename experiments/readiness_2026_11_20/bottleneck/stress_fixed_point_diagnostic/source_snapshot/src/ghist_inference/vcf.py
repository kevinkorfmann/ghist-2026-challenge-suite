"""VCF checks shared by fixed-sample demographic spectrum extractors."""


def validate_position(record, contig_lengths):
    length = contig_lengths.get(record.contig)
    if (
        record.contig not in contig_lengths
        or record.pos < 1
        or (length is not None and record.pos > length)
    ):
        raise ValueError(f"Invalid 1-based VCF coordinate {record.contig}:{record.pos}")


def diploid_alleles(record, samples, *, require_complete=False):
    alleles = []
    for sample in samples:
        genotype = record.samples[sample].get("GT")
        if genotype is None:
            genotype = (None, None)
        if len(genotype) != 2 or any(a not in (None, 0, 1) for a in genotype):
            raise ValueError(
                f"Expected biallelic diploid GT for {sample} at {record.contig}:{record.pos}"
            )
        if require_complete and None in genotype:
            raise ValueError(
                "Missing genotype in fixed-sample spectrum; use a projected spectrum and callable exposure instead"
            )
        alleles.extend(a for a in genotype if a is not None)
    return alleles
