stringtie ../mapping/SL1.bam  -e -B -A sl1.gene.abundance -G ../genome/annotation.all_transcripts.all_features.peq.gff3 -p 20 -o sl1/sl1.gtf
stringtie ../mapping/SL2.bam  -e -B -A sl2.gene.abundance -G ../genome/annotation.all_transcripts.all_features.peq.gff3 -p 20 -o sl2/sl2.gtf
stringtie ../mapping/SL3.bam  -e -B -A sl3.gene.abundance -G ../genome/annotation.all_transcripts.all_features.peq.gff3 -p 20 -o sl3/sl3.gtf
stringtie ../mapping/SF1.bam  -e -B -A sf1.gene.abundance -G ../genome/annotation.all_transcripts.all_features.peq.gff3 -p 20 -o sf1/sfl.gtf
stringtie ../mapping/SF2.bam  -e -B -A sf2.gene.abundance -G ../genome/annotation.all_transcripts.all_features.peq.gff3 -p 20 -o sf2/sf2.gtf
stringtie ../mapping/SF3.bam  -e -B -A sf3.gene.abundance -G ../genome/annotation.all_transcripts.all_features.peq.gff3 -p 20 -o sf3/sf3.gtf
