ls CLEANREADS/*_1.clean.fq.gz | parallel -j4 'echo {/} $(zcat {} | awk "{s++}END{print s/4}")' > readCounts.tsv
