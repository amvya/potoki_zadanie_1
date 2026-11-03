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
# WARNING: initially runs the pipeline up to the rule merge_flt
rule all:
    input:
        'output/tree.newick',
        'output/tree.png'

# Rule blastn: for each genome perform a nucletide BLAST search using:
# - a single-sequence (a fragment of gmk gene) FASTA query file
# - a multiple-sequence FASTA file (one staphylococcal genome assembly contigs)
# Uncompress the gzipped genome assembly file and redirect the output
# to blastn's stdin while asking blastn to read the subject FASTA sequence
# from that stream ('-' instead of a subject file path)
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
# rule merge_mlst:
    #TODO: implement the rule
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
# rule unique:
    #TODO: implement the rule
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
# rule tree:
    #TODO: implement the rule
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

