python3 main-script-targeted-caller.py \
    -i "example_input/sample_list.txt" \
    -b "example_input/control_list.txt" \
    -p "example_input/parsed_lane_output" \
    -wm targeted-caller/NCCN_SNV_indel.txt \
    -wf targeted-caller/NCCN_fusions.txt \
    -o example_output\
    -s 1234 \
    -g /silo5/cfRNA/indexes_v3.3/hg19_ERCC.fa \
    -bt bedtools\
    $@  