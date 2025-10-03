GVCFS=/data3/projects/mirage/as-internship-files-26022021/hb3Mdh2/gvcfs/*.g.vcf
for f in $GVCFS
  do
   bgzip $f
   tabix -p vcf $f.gz
 done

