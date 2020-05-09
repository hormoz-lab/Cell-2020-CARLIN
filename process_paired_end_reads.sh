#!/bin/bash
#
# Uses PEAR to created and compressed PE reads from
# from input and output locations specified in SAMPLE_SHEET 

display_usage() { 
   
    echo "This script uses PEAR to merge single-end reads from RAW_PATH"\
    "and places the zipped results in PROCESSED_PATH. It does this for"\
    "all files in the subdirectory structure provided in SAMPLE_SHEET."\
    "If -SRA is optionally specified, the Illumina subdirectory structure"\
    "is ignored and a flattened subdirectory structure obtained from"\
    "SRA-tools fastq-dump is used instead."    
    
    echo -e "\nUsage:\n\n$0 PEAR_PATH/bin/pear SAMPLE_SHEET RAW_PATH PROCESSED_PATH [-SRA]\n"
}
	
if [  $# -lt 4 ]; then 
    display_usage
    exit 1
fi

SRA=0

if [ $# == 5 ]; then
 
    if [ $5 != "-SRA" ]; then
        display_usage
        exit 1
    fi

    SRA=1

fi 


if [[ ( $# == "--help") ||  $# == "-h" ]]; then 
    display_usage
    exit 0
fi


PEAR_APP=$1
SAMPLE_SHEET=$2
RAW_PATH=$3
PROCESSED_PATH=$4

mkdir -p $PROCESSED_PATH


cat $SAMPLE_SHEET | while read line; do
    
    line=( ${line//,/ } )
  
    if [ "${line[6]}" = "PEAR" ]; then
        echo "${line[7]}" 
        mkdir -p $PROCESSED_PATH/${line[7]}

        if [[ $SRA = 1 ]]; then
            F_FILE=$RAW_PATH/${line[10]}_1.fastq 
            R_FILE=$RAW_PATH/${line[10]}_2.fastq
            gunzip -k $F_FILE.gz
            gunzip -k $R_FILE.gz
        else
            F_FILE=$RAW_PATH/${line[2]}/${line[5]}_R1_001.fastq.gz 
            R_FILE=$RAW_PATH/${line[2]}/${line[5]}_R2_001.fastq.gz 
        fi

        echo $F_FILE
        echo $R_FILE

        $PEAR_APP -f $F_FILE \
                  -r $R_FILE \
                  -o $PROCESSED_PATH/${line[7]}/PE \
                  --min-overlap 1 \
                  --min-assembly-length 1 \
                  --min-trim-length 1 \
                  --quality-threshold 30 \
                  --threads 4 \
                  > $PROCESSED_PATH/${line[7]}/pear.log

        cd $PROCESSED_PATH/${line[7]}
        gzip *.fastq        
        cd $OLDPWD

        if [[ $SRA = 1 ]]; then
            rm $F_FILE
            rm $R_FILE
        fi

    fi

done

