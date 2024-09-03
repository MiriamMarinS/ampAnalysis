name_dir=$1

for file in ${name_dir}*_R1*.fastq
do
    R1=$(basename $file)
    echo 'Processing sample: '${R1}
    R2=$(echo $R1 | sed 's/R1/R2/g')
    fileR2=$(echo $file | sed 's/R1/R2/g')
    barcode=$(more ./barcodes.txt | grep $R1 | cut -d$'\t' -f 3)
    cutadapt -g $barcode -G $barcode -o ./$R1 -p ./$R2 $file $fileR2
done