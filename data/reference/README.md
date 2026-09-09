# Reference Files

## GRCh38 (hg38) to GRCh37 (hg19) LiftOver Chain File

The ovarian cancer GWAS summary statistics used in this study are
provided in the GRCh38 (hg38) genome build. Before performing the
pleiotropy analysis, the ovarian cancer variants are converted from
GRCh38 (hg38) to GRCh37 (hg19) using the UCSC LiftOver chain file:

`hg38ToHg19.over.chain`

### Download

Download the `hg38ToHg19.over.chain` file from the UCSC Genome Browser
LiftOver resources.

After downloading, place the file in this directory:

```text
data/reference/hg38ToHg19.over.chain
