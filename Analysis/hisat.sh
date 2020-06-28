for i in `cat ../sample`;do
echo $i;
hisat2 -p 10 --dta --rg-id $i --rg SM:$i --rna-strandness RF --fr -x ../genome/pha -1 ../trimmed/${i}_1P.fq.gz -2 ../trimmed/${i}_2P.fq.gz -S ${i}.sam
samtools view -bS ${i}.sam > ${i}.bam;samtools sort -@ 10 ${i}.bam ${i}.sorted;samtools index ${i}.sorted.bam;rm ${i}.sam
done
