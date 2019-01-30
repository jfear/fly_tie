from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from snakemake import shell

# Objects passed in by snakemake
srx = snakemake.wildcards.sample
sampletable = pd.read_csv(snakemake.config['sampletable'], sep='\t', index_col=0)
oname1 = snakemake.output.r1
oname2 = snakemake.output.r2

# Download in a tmpdir
TMPDIR = TemporaryDirectory()


def main():
    srrs = get_srrs()
    fnames = [fastq_dump(srr) for srr in srrs]
    r1 = [read[0] for read in fnames]
    r2 = [read[1] for read in fnames]
    cat_files(r1, oname1)
    cat_files(r2, oname2)


def get_srrs():
    srrs = sampletable.loc[srx, 'Run']

    if isinstance(srrs, str):
        return [srrs]

    return srrs.tolist()


def fastq_dump(srr):
    shell(
        'cd {TMPDIR.name} '
        '&& fastq-dump '
        '{srr} '
        '--gzip '
        '--split-files '
        '&& touch {srr}_1.fastq.gz '
        '&& touch {srr}_2.fastq.gz '
    )

    return f'{srr}_1.fastq.gz', f'{srr}_2.fastq.gz'


def cat_files(fnames, oname):
    # remove file if it is there
    if Path(oname).exists():
        shell('rm {oname}')

    # Iterate over FASTQs and concat
    for fname in fnames:
        shell('cat {TMPDIR.name}/{fname} >> {oname}')


if __name__ == '__main__':
    main()
