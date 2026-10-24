# yruler

Measure the Y chromosome satellite arrays from sequencing reads.

yruler counts marker k-mers built from complete Y assemblies and reports:

- length of the Yq12 array and of its DYZ1 and DYZ2 parts;
- length and composition of the DYZ3 centromeric array;
- Y dosage, the fraction of cells that carry the Y;
- quality flags for low Y dosage and loss of part of the Y.

It reads aligned (BAM/CRAM) or unaligned genomes (FASTQ). For aligned genomes, it fetches only the reads that can carry markers.

