#!/bin/bash -ue
cp multiqc/* .
echo "custom_logo: $PWD/logo.png" >> multiqc_config.yaml
multiqc -v .
