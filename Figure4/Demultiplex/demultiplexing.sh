#!/bin/bash

#######################################
#     Get all needed variables        #
#######################################

# xcalibr-master directory path
#------------------------------
echo -n "

    XCALIBR DIRECTORY 
    -----------------
    
What is your tool path (place where you git clone xcalibr) ?

Enter the xcalibr-master directory path. 
    Ex: /Users/hadjabed/Documents/Curie/TPerie/xcalibr-master 

";
read;
XCAL_TOOLPATH="${REPLY}"

while [ ! -d "$XCAL_TOOLPATH" ]
do
    echo -n "
    This directory does not exists.
    Try again:

";
read;
XCAL_TOOLPATH="${REPLY}"
done

# Fastq directory 
#----------------
echo -n "

    FASTQ DIRECTORY 
    ---------------
    
What is your fastq directory path ?

    Ex: /Users/hadjabed/Dropbox (Team_Perie)/Louisa_Hadj/Fastq

";
read;
IN_DIR="${REPLY}"

while [ ! -d "$IN_DIR" ]
do
    echo -n "
    This directory does not exists.
    Try again: 
    
";
read;
IN_DIR="${REPLY}"
done


# LIBRARY
#--------
echo -n "

    LIBRARY INFORMATION 
    -------------------
    
What is the library used ?

If it is the old one type: old
If it is the new one type: new

NB: write your answer extactly the same way than propositions. 

";
read;
rep_refLib="${REPLY}"

# tant que vrai
while [[ "$rep_refLib" != "old" ]] && [[ "$rep_refLib" != "new" ]]
do 
    echo -n "
    You did not type old or new.

    Again, if you use the old library type: old
    If you use the new one type: new

    ";
    read;
    rep_refLib="${REPLY}"
done


if [[ $rep_refLib == "old" ]]
then 
    # OLD LIB
    LIB_DIR="./inputs/Libraries/old_lib"
    REF_LIB="reference_barcode_library_lentirus.txt"
    INDEX_REF="./inputs/Libraries/old_lib/index_library_384.txt"
    CONST_POS="9-24"
    #ROW="3"
    DESIGN="X8N16Y15"
else
    # NEW LIB
    INDEX_REF="./inputs/Libraries/new_lib/plate_index_library.txt"
    CONST_POS="9-27"
    ROW="3,6"

    echo -n "
    What is the barcode length ?

    If it is 20bp type: 20
    If it is 21bp type: 21
    
    NB: write your answer extactly the same way than propositions. 

    ";
    read;
    rep_size="${REPLY}"

    # tant que vrai
    while [[ "$rep_size" != 20 ]] && [[ "$rep_size" != 21 ]]
    do 
        echo -n "
        You did not type 20 or 21.

        Again, if your barcode is 20 base long type: 20
        Or if your barcode is 21 base long type: 21

        ";
        read;
        rep_size="${REPLY}"
    done

    if [[ "$rep_size" == 20 ]]
    then
        DESIGN="X8N19Y20"
        LIB_DIR="./inputs/Libraries/new_lib/Ref_Lib_20bp"
    else
        DESIGN="X8N19Y21"
        LIB_DIR="./inputs/Libraries/new_lib/Ref_Lib_21bp"
    fi

    echo -n "
    Which library did you used ?

    If it is LG21 type: LG21
    If it is LG22 type: LG22
    If it is LT21 type: LT21
    If it is LT22 type: LT22
    
    NB: write your answer extactly the same way than propositions. 

    ";
    read;
    rep_lib="${REPLY}"

    # tant que vrai
    while [[ "$rep_lib" != "LG21" ]] && [[ "$rep_lib" != "LG22" ]] && [[ "$rep_lib" != "LT21" ]] && [[ "$rep_lib" != "LT22" ]]
    do 
        echo -n "
        You did not type LG21, LG22, LT21 or LT22.

        Again, which of the above propositions correspond to your library ?

        ";
        read;
        rep_lib="${REPLY}"
    done

    REF_LIB=$rep_lib"_filtered.fa"
fi

SAMPLE_DESCRIPTION=$(ls ./inputs/sampleDescription)
SAMPLE_DESCRIPTION_PATH="./inputs/sampleDescription/"$SAMPLE_DESCRIPTION

OUT_DIR_xcalibr="./outputs/xcalibr"
OUT_DIR_final="./outputs/final_matrices"
TMP_DIR="./outputs/tmp/"

###########################
#     1. xcalibr          #
##########################


## Add xcalibr app location to your PATH 
PATH=$PATH:"$(echo $XCAL_TOOLPATH)"
PATH=$PATH:"$(echo $XCAL_TOOLPATH)"/bin
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

for fastq in `ls "$(echo $IN_DIR)"` 
do
    echo $fastq
    prefix=`basename $fastq .R1.fastq.gz`
    mkdir "$OUT_DIR_xcalibr/$prefix"
    FILE_RES=`echo "$OUT_DIR_xcalibr"/"$prefix" `
    ## 1) 'hash' create a binary and counted format.
    gzip -cd "$(echo $IN_DIR)"/$fastq | cat | xcalibr hash "$FILE_RES"/"$prefix".R1.bin
    # 2) locate and extract the barcodes
    xcalibr extract --template "$DESIGN" --matchX "$INDEX_REF" --matchY "$LIB_DIR"/$REF_LIB --printemptyX "$FILE_RES"/"$prefix".R1.bin > "$FILE_RES"/"$prefix"-matrix.txt
    ## 3) 'analyze' read the hashed representation of the sequences and analyze the given position. Print top occurring sequences in that range.
    # OLD LIB: constant part (botlib) = 9-24 
    # NEW LIB: constant part (toplib) = 9-28
    xcalibr analyze "$FILE_RES"/"$prefix".R1.bin $CONST_POS > "$FILE_RES"/"$prefix"_constantPart-analysis.txt
    # only new lib has plate index:
    if [[ $rep_refLib == "new" ]]
    then 
        xcalibr analyze "$FILE_RES"/"$prefix".R1.bin 1-8 > "$FILE_RES"/"$prefix"_index-analysis.txt
    fi
done


###########################################
#     2. Reformat xcalibr outputs         #
###########################################

### 2.1 - Create final matrix ###
#-------------------------------#

analysisID=$(basename "$SAMPLE_DESCRIPTION_PATH" ".sampleDescription.txt")
echo $analysisID > analysis_ID

for matrix in $OUT_DIR_xcalibr/*/*-matrix.txt
do
    fileName=$(basename $matrix -matrix.txt)
    fileNum=$(echo $fileName | sed "s/$analysisID//")
   
    # NEW
    if [[ $rep_refLib == "new" ]]
    then 
        sample_name=`grep $fileNum\| "$SAMPLE_DESCRIPTION_PATH"  | cut -d '|' -f2 `
        cut -f1-5 $matrix > $TMP_DIR/$sample_name.csv
    fi
    # OLD
    if [[ $rep_refLib == "old" ]]
    then 
        for line in $(cat "$SAMPLE_DESCRIPTION_PATH")
        do
            sampleIndex=$(echo $line | cut -d "|" -f1 )
            sample_name=$(echo $line | cut -d "|" -f2 )
            sed "s/$sampleIndex/$sample_name/" $matrix > tmp
            mv tmp $matrix
        done
        cp $matrix $TMP_DIR/
    fi


done

Rscript bin/create_finalMatrix.R $analysisID $rep_refLib


### 2.2 Create constant part & index summary ###
#----------------------------------------------#
i=1
for file in $OUT_DIR_xcalibr/*/*-analysis.txt
do
    # if file name contains "constantPart"
    if [[ "$file" =~ constantPart ]]
    then
        # Get colnames = sequence of the constant part
        if [ "$i" -eq 1 ]
        then 
            seq=$(sed -n 3p $file | awk '{print $1}')
            echo "%" $seq  > results_constantPart
            echo ' ' > names
            i=2
        fi

        # Get constant part percentages
        percent=$(sed -n 3p $file | awk '{print $3}')
        echo $percent | sed $'s/ /\t/g' >> results_constantPart

        # Get row names of constant parts results. 
        # For the new lib, each line is a sample name, for the old lib each line is an experiement. 
        fileName=$(basename $file _constantPart-analysis.txt)
        if [[ $rep_refLib == "new" ]]
        then 
            fileNum=$(echo $fileName | sed "s/$analysisID//")
            sample_name=`grep $fileNum\| "$SAMPLE_DESCRIPTION_PATH" | cut -d '|' -f2`
            echo $sample_name >> names
        else
            echo $fileName >> names
        fi

    fi
    # if file name contains "index" (only new lib)
    if [[ "$file" =~ index ]]
    then
        # Get index part sequences and percentages
        if [ "$i" -eq 2 ]
        then 
            seq=$(sed -n "$ROW"p $file | awk '{print $1}')
            echo '%'$seq | sed $'s/ /\t/g' > results_index
        fi
        i=3
        percent=$(sed -n "$ROW"p $file | awk '{print $3}')
        echo $percent | sed $'s/ /\t/g' >> results_index
    fi

done

## Constant part
paste names results_constantPart > $OUT_DIR_final/constantPart_summary.tsv

## Plate index part
if [[ $rep_refLib == "new" ]]
then 
    paste names results_index > $OUT_DIR_final/index_summary.tsv
    rm results_index analysis_ID
fi

### Clean folder
rm names results_constantPart
