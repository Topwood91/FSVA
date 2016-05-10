# FSVA
README for FSVA

Installation
--------------
Download all these files(fsva.cpp, index.cpp, ssw.h, ssw.c, ssw_cpp.h, ssw_cpp.cpp, radix_sort.cpp, makefile) and put them in the same folder.
Then, typing
make
to build. If the build is successful, two new executable files are created:index and fsva.

Index
--------------
To do alignment, index is needed first. The command of index is

./index [options] <reference.fa>

options:
        -r INT length of read. Using this argument, length of seed
               can be calculated by program automatically. [150]
        -l INT length of seed. [31]

As FSVA will choose different seed length to fit different read length, we need know the read length when indexing. Users can set the seed length by themselves, but it is not recommended. Notice that don't use -r and -l at the same time.

Alignment
--------------
After indexing, users can use fsva to do alignment. The command of alignment is

./fsva [options] <reference.fa> <read1.fastq> [read2.fastq]

options:
        -t INT number of threads. [1]
        -u INT threshold of unrepresentative seed. [450]
        -l INT length of seed. If you use this argument when making
               index, you should use this argument here, and their
               value should be equal. [31]
        -r INT length of read. If you use this argument when making
               index, you should use this argument here, and their
               value should be equal. [150]

If you do alignment on single-end data, your command may like this:

./fsva [options] hg19.fa read1.fq > out.sam

and if you do alignment on pair-end data, your command may like this:

./fsva [options] hg19.fa read1.fq read2.fq > out.sam

Here, -u represents the threshold of unrepresentative seed. If you use a small value, more seeds will be droped. As a result, the speed of program will increase while accuracy will decline. 
To guarantee the program work correctly, if you use -r or -l when index, you should use it here, and there value should be equal.
