#! /bin/bash
##################################################
##Date_of_creation=9/10/2020
##Tested successfully on 13/10/2020

##############CREATING DIRECTORY VARIABLES########
function producing_a_vcf_from_bam() {

################GATK_function######################
#Input bam file to be exported from the master file
local INPUTFILE="$1"

#########################################################
###################DATA-Preprocessing####################

#Creating Sequence Dictionaries
java -jar /usr/local/picard-tools-2.5.0/picard.jar CreateSequenceDictionary \
     R="$REFERENCE" \
     O=$INPUTFILE.dict
echo "SequenceDictionary created for sample $INPUTFILE"


#complete the header of each bam file using the AddOrReplaceReadGroups from Picard
file_name="${INPUTFILE##*/}"
onlyFileName="${file_name%.*}" #remove.bam for SM information
onlyFileName="${onlyFileName%.*}" #remove .sorted SM information
java "-Xmx1G" "-Xms1G" "-XX:ConcGCThreads=1" -jar /usr/local/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups \
        I=$INPUTFILE \
        O=$INPUTFILE.fix \
        LB=WGA \
        PL=illumina \
        PU=NA \
        SM=$onlyFileName
local fate=$?
if [[ $fate -ne 0 ]]
        then
		return  [1]
		echo "$INPUTFILE encountered a problem while AddOrReplaceReadGroups Processing"
else
	echo "AddOrReplaceReadGroups PROCESSED $f"
fi

#Renaming the output file so that it acts as input for the next step
mv $INPUTFILE.fix $INPUTFILE

#index each fixed bam file using samtools
samtools index $INPUTFILE
echo "index PROCESSED $f"

#Using the GATK haplotype Caller
gatk --java-options "-Xmx4G"  HaplotypeCaller \
	-R "$REFERENCE" \
	-I $INPUTFILE \
	--emit-ref-confidence GVCF \
	--pcr-indel-model NONE \
	--sample-ploidy 1 \
	--max-alternate-alleles 2 \
	--output $INPUTFILE.g.vcf

#Moving the gvcfs to a subfolder
mv "$INPUTFILE".g.vcf "$OUTPUTDIR"
echo "GVCF produced for $INPUTFILE in $OUTPUTDIR"
}

"$@"
