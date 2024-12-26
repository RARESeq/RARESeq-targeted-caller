python3 parse-lane-outputs.py \
    -i "example_input/sample_list.txt" \
    -b "example_input/control_list.txt" \
    -l /silo5/cfRNA/v3.3 \
    -wm targeted-caller/NCCN_SNV_indel.txt \
    -wf targeted-caller/NCCN_fusions.txt \
    -o example_input\
    -s 1234 \
    -g /silo5/cfRNA/indexes_v3.3/hg19_ERCC.fa \
    -bt bedtools\
    $@  