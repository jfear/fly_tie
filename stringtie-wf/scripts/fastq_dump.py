from pathlib import Path
from tempfile import TemporaryDirectory

from snakemake import shell

# Objects passed in by snakemake
sampletable = snakemake.params.sampletable
srx = snakemake.wildcards.sample
oname = snakemake.output[0]

# Download in a tmpdir
TMPDIR = TemporaryDirectory()


def main():
    srrs = get_srrs()
    layout = get_layout()
    fnames = [fastq_dump(srr) for srr in srrs]
    r1 = [read[0] for read in fnames]
    r2 = [read[1] for read in fnames]

    cat_files(r1, oname)

    if layout == 'PE':
        cat_files(r2, oname.replace('R1', 'R2'))


def get_srrs():
    srrs = sampletable.loc[srx, 'Run']

    if isinstance(srrs, str):
        return [srrs]

    return srrs.tolist()


def get_layout():
    layouts = sampletable.loc[srx, 'layout']
    if isinstance(layouts, str):
        return layouts
    return layouts.tolist()[0]


def fastq_dump(srr):
    shell(
        'cd {TMPDIR.name} '
        '&& fastq-dump '
        '{srr} '
        '--gzip '
        '--split-files'
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
