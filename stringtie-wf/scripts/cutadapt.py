from snakemake import shell

fastqs = snakemake.input
output = snakemake.output
extra = snakemake.params.extra
log = snakemake.log


def main():
    paried = len(fastqs) == 2
    oname1 = output[0]
    oname2 = oname1.replace('R1', 'R2')

    if paried:
        shell(
            "cutadapt "
            "{extra} "
            "{fastqs[0]} "
            "{fastqs[1]} "
            "-o {oname1} "
            "-p {oname2} "
            "&> {log}"
        )
    else:
        shell(
            "cutadapt "
            "{extra} "
            "{fastqs[0]} "
            "-o {oname1} "
            "&> {log}"
        )


if __name__ == '__main__':
    main()
