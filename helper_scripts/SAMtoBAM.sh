
# Use samtools in shell script create bam files
/Users/caryodonnell/unixtools/bin/samtools view -S -b $1/$2.sam > $1/$2.bam
/Users/caryodonnell/unixtools/bin/samtools sort $1/$2.bam -o $1/$2.sorted.bam