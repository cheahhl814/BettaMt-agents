Usage: circlator <command> [options] <required arguments>

To get minimal usage for a command use:
circlator command

To get full help for a command use one of:
circlator command -h
circlator command --help


Available commands:

all        Run mapreads, bam2reads, assemble, merge, clean, fixstart
mapreads   Map reads to assembly
bam2reads  Make reads from mapping to be reassembled
assemble   Run assembly using reads from bam2reads
merge      Merge original assembly and new assembly made by assemble
clean      Remove small and completely contained contigs from assembly
fixstart   Change start position of circular sequences
minimus2   Run the minimus2 based circularisation pipeline
get_dnaa   Download file of dnaA (or other of user's choice) genes
progcheck  Checks dependencies are installed
test       Run Circlator on a small test set
version    Print version and exit
