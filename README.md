## Important Notes
\* GFR means GeneFuseRust
#### - Caveats in fastq reading
~~GFR reads lines from fastq file, with maximum read bytes 1000.   
If a line has more than 1000 bytes, the result files may have critical errors and the program does not warn or is terminated.~~  
-> 240120 fixed.
