# =============================================================================
# MLST-based phylogeny of staphylococcal genomes (Snakemake workflow)
#
# Overview of the workflow (one chain of rules per genome assembly):
#
#   blastn -> filter -> merge_flt -> mlst -> merge_mlst -> unique -> tree
#
#   1. blastn     search every MLST locus query against every genome assembly
#   2. filter     keep only full-length, gapless hits (complete allele copies)
#   3. merge_flt  merge the filtered hits of all loci into one table per genome
#   4. mlst       assign allele numbers and a sequence type (ST) per genome
#   5. merge_mlst merge per-genome MLST results into one table
#   6. unique     extract one concatenated locus sequence per ST (FASTA)
#   7. tree       build a neighbour-joining (NJ) tree and draw it
#
# Wildcards used throughout the workflow:
#   {seqid} - name of a query (locus) file, input/queries/{seqid}.fna
#   {asmid} - genome assembly accession, <genomesdir>/{asmid}.fna.gz
#
# Usage: see README.md (snakemake --use-conda --cores <N>)
# =============================================================================


# Indicate the config file
configfile:
    'config.yml'

# Extract query sequence names (sequence IDs) from file paths
# to use them as IDs for query files
seqids, = glob_wildcards('input/queries/{seqid}.fna')

# Extract genome assembly accession numbers (assembly IDs) from file paths
# to use them as IDs for genome files
asmids, = glob_wildcards(config['genomesdir'] + '/{asmid}.fna.gz')

# Rule all: collect the final output (multiple files)
# The target files are the NJ tree in Newick format and its image; Snakemake
# derives and runs all the preceding rules needed to produce them
rule all:
    input:
        'output/tree.newick',
        'output/tree.png'

# Rule blastn: for each genome perform a nucleotide BLAST search using:
# - a single-sequence (a fragment of an MLST locus, e.g. the gmk gene)
#   FASTA query file
# - a multiple-sequence FASTA file (one staphylococcal genome assembly contigs)
# Uncompress the gzipped genome assembly file and redirect the output
# to blastn's stdin while asking blastn to read the subject FASTA sequence
# from that stream ('-' instead of a subject file path)
# The output is a tab-separated table (-outfmt 6) with the columns listed
# in params.cols; sseq (aligned subject sequence) is later used as the allele
rule blastn:
    conda:
        'envs/powb-blast.yml'
    params:
        evalue = 0.01,
        cols   = 'qseqid qlen sseqid length qcovs pident evalue sseq'
    input:
        query = 'input/queries/{seqid}.fna',
        subj  = config['genomesdir'] + '/{asmid}.fna.gz'
    output:
        'output/blastn/blastn_{asmid}_{seqid}.tsv'
    log:
        'log/blastn/blastn_{asmid}_{seqid}.log'
    shell:
        '''gunzip -c {input.subj} 2> {log} |      \
           blastn -evalue         {params.evalue} \
                  -outfmt      "6 {params.cols}"  \
                  -query          {input.query}   \
                  -subject        -               \
                  -out            {output}        \
                   >> {log} 2>&1
        '''

# Rule filter: for each search retain 100% query coverage hits (alignments)
# of the same length as the query sequence (gapless alignment)
# (--qcovs sets the minimum query coverage, --eqlen requires the alignment
# length to equal the query length; both are handled in scripts/filter.py)
rule filter:
    conda:
        'envs/powb-python.yml'
    params:
        qcovs = 100,
        cols  = rules.blastn.params.cols
    input:
        rules.blastn.output
    output:
        'output/filter/filtered_{asmid}_{seqid}.tsv'
    log:
        'log/filter/filter_{asmid}_{seqid}.log'
    shell:
        '''scripts/filter.py --eqlen                    \
                             --asmid  {wildcards.asmid} \
                             --qcovs  {params.qcovs}    \
                             --cols  "{params.cols}"    \
                             --input  {input}           \
                             --output {output}          \
                               > {log} 2>&1
        '''

# Rule merge_flt: join all filtered BLAST result files into one using
# a short Bash script (a list of file path on input)
# One merged table is created per genome (all loci of that genome); the header
# line is taken from the first file and the header lines of the others are skipped
rule merge_flt:
    params:
        mask = rules.filter.output[0].replace('{seqid}', '*')
    input:
        expand(rules.filter.output, seqid=seqids, allow_missing=True)
    output:
        'output/merge_flt/merged_{asmid}.tsv'
    log:
        'log/merge_flt/merge_{asmid}.log'
    shell:
        '''files=({params.mask})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail -n +2 ${{file}}
           done >> {output} 2>> {log}
        '''

# Rule mlst: perform MLST analysis, one genome at a time
# Using the merged hits of a genome, the alleles are matched against
# the allele database (input/alleles.fna) and the allele profile is looked up
# in the MLST profile table (input/bigsdb.tsv) to assign the sequence type (ST)
# See scripts/mlst.py for details
rule mlst:
    conda:
        'envs/powb-python.yml'
    input:
        alleles = expand(rules.merge_flt.output, seqid=seqids, allow_missing=True), 
        allelesdb = 'input/alleles.fna',
        bigsdb = 'input/bigsdb.tsv'
    output:
        'output/mlts/mlst_{asmid}.tsv'
    log:
        'log/mlst/mlst_{asmid}.log'
    shell:
        '''
        scripts/mlst.py \
            --asmid {wildcards.asmid} \
            --allelesdb {input.allelesdb} \
            --mlstdb {input.bigsdb} \
            --alleles {input.alleles} \
            --output {output} \
            > {log} 2>&1
        '''

# Rule merge_mlst: join all MLST results files into one using a short Bash
# script (a list of file path on input)
# The result is a single table with one row per genome (the header line is
# taken from the first file). NOTE: despite the .sh extension of the output
# file name, the file contains a tab-separated table, not a script
rule merge_mlst:
    input: expand(rules.mlst.output, asmid=asmids, allow_missing=True), 
    output:
        'output/merge_mlst.sh'
    log:
        'log/merge_mlst.log'
    shell:
        '''files=({input})
           head -1 ${{files[0]}} > {output} 2> {log}
           for file in "${{files[@]}}"; do
               tail +2 $file >> {output} 2>> {log}
           done
        '''

# Rule unique: from merged BLAST results extract unique merged MLST loci
# sequences (i.e., one sequence per ST) to a FASTA file
# Used as the sequence set for tree construction (see scripts/unique.py)
rule unique:
    conda: 'envs/powb-python.yml'
    input: rules.merge_mlst.output
    output:
        'output/unique.fna'
    log:
        'log/unique.log'
    shell:
        '''
        scripts/unique.py --input {input} --output {output}
        '''

# Rule tree: take the merged MLST results and extracted merged MLST loci
# sequences, and then compute and visualise a phylogenetic tree (NJ tree)
# Outputs: the tree in Newick format and its graphical representation (PNG)
# (see scripts/tree.py)
rule tree:
    conda: 'envs/powb-python.yml'
    input: 
        inputunique = rules.unique.output,
        inputmerge = rules.merge_mlst.output
    output:
        outtree = 'output/tree.newick',
        outvis = 'output/tree.png'
    log:
        'log/tree.log'
    shell:
        '''
        scripts/tree.py \
            --mlst {input.inputmerge} \
            --seqs {input.inputunique} \
            --outtree {output.outtree} \
            --outvis {output.outvis} \
            > {log} 2>&1 
        '''

