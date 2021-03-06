#!/bin/bash
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
   echo "o     output.depth"
   echo
   echo "example: GGmmDepth.sh -i input.bam -o output.depth"
}

############################################################
while getopts ":hio:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
         i) # Enter input
         input=$OPTARG;;
         o) # Enter output
         out=$OPTARG;;
     \?) # Invalid option
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
samtools view -H input > Xheader.txt
samtools view -h input | grep -v "NM:i:0" | grep -E "MD:Z:.*([0-9]{1,3}[CG][0-9]{1,3})" | cat Xheader.txt - |  samtools view -Sb | samtools depth -aa - > out
rm Xheader.txt
############################################################
