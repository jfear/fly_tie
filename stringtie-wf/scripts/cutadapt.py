import pandas as pd
from snakemake import shell

fastq1 = snakemake.input.r1
fastq2 = snakemake.input.r2

oname1 = snakemake.output.r1
oname2 = snakemake.output.r2

sampletable = pd.read_csv(snakemake.config['sampletable'], sep='\t', index_col=0)

extra = snakemake.params.extra
log = snakemake.log


def main():
    layout = get_layout()

    if layout == 'PE':
        shell(
            "cutadapt "
            "{extra} "
            "{fastq1} "
            "{fastq2} "
            "-o {oname1} "
            "-p {oname2} "
            "&> {log}"
        )
    else:
        shell(
            "cutadapt "
            "{extra} "
            "{fastq1} "
            "-o {oname1} "
            "&> {log} "
            "&& touch {oname2}"
        )


def get_layout():
    srx = snakemake.wildcards.sample
    layout = sampletable.loc[srx, 'layout']
    if isinstance(layout, str):
        return layout
    return layout.tolist()[0]


if __name__ == '__main__':
    main()
