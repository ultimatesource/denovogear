#!/bin/bash
npm run build

python2.7 tools/build.py --rscript_location /usr/bin/Rscript --ped_file test/data/test.ped --source_dir dist --dng_output_file_path test/data/dng_test_output.vcf --output_file_path mutmap.html
