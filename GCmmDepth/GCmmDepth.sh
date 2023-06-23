#!/bin/bash
#### function modified on 06/22/2022 
#### R.S.
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This function extracts the depth of coverage across the genome after filtering out sequencing reads without at least one mismatch at any G:C position"
   echo "! requires samtools !"
   echo 
   echo "Syntax: GGmmDepth.sh [h] -i input -o output" 
   echo "options:"
   echo "h     Print this Help."
   echo "i     input.bam"
   echo "o     output.depth; use "-" for stdout"
   echo
   echo "example: GGmmDepth.sh -i input.bam -o output.depth"
}

############################################################
while getopts ":i:o:h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # input file
         input=$OPTARG;;
      o) # output file
         output=$OPTARG;;
     \?) # incorrect option
         echo "Error: Invalid option"
         exit;;
   esac
done
############################################################
# Main Function                                            #
############################################################
if [[ $# -eq 0 ]] ; then
    echo 'Error: No argument supplied'
    echo 'use -h option for help'
    exit 1
fi

## test if input file exists
if [ -f $input ]; then
   echo "Getting depth of coverage from reads with at least one mismatch at any G:C position from $input"
else
    echo "Error: Input file not found"
    exit 1
fi

### test if $output argument is missing
if [ -z "$output" ]; then
    echo "Error: Output file not specified"
    exit 1
fi
## test is $output is equal to "-"
if [ $output == "-" ]; then
samtools view -H $input > Xheader.txt
samtools view -h $input | grep -v "NM:i:0" | grep -E "MD:Z:.*([0-9]{1,3}[CG][0-9]{1,3})" | cat Xheader.txt - |  samtools view -Sb | samtools depth -aa - 
rm Xheader.txt
fi


if [ -f $output ]; then
   echo "Output $output already exists"
   echo "Will overwrite $output"
else
   echo "Depth of coverage will be written to $output"
fi

samtools view -H $input > Xheader.txt
samtools view -h $input | grep -v "NM:i:0" | grep -E "MD:Z:.*([0-9]{1,3}[CG][0-9]{1,3})" | cat Xheader.txt - |  samtools view -Sb | samtools depth -aa - > $output
rm Xheader.txt

#### test if output file is created
if [ -f $output ]; then
    echo "Completed; Depth of coverage written to $output"
else
    echo "Error: Depth of coverage file not created"
    exit 1
fi
############################################################
