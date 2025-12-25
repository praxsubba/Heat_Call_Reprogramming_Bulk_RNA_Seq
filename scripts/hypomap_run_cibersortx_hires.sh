apptainer exec \
  -B /path/to/input:/src/data \
  -B /path/to/output/Hypomap_output:/src/outdir \
  /path/to/hires_latest.sif \
  /src/CIBERSORTxHiRes \
    --username YOUR_EMAIL_HERE \
    --token YOUR_TOKEN_HERE \
    --mixture /src/data/CIBERSORTx_Mixtures_Adjusted.txt \
    --sigmatrix /src/data/Hypomap_custom_signature_matrix.txt \
    --QN FALSE
