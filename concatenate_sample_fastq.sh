for s in $(cat samples_sequenced_thrice.txt);
do
    echo $s 
    
    cat data/fastq/from_ora/${s}.1.fastq.gz >> \
    data/fastq/${s}.1.fastq.gz 

    cat data/fastq/from_ora/${s}.2.fastq.gz >> \
    data/fastq/${s}.2.fastq.gz
done

