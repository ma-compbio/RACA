RACA (Reference-Assisted Chromosome Assembly)
=============================================

This is the version used in the original [PNAS][PNAS] paper.


How to compile?
---------------

Just type `make` to compile the RACA package.


How to run?
-----------

### Configuration file

RACA requires a single configuration file as a parameter. 
The configuration file has all parameters that are needed for RACA.

Example dataset below has a sample configuration file 'params.txt' 
and all other parameter files which are self-explanatory.
(also refer to the data directory.) 

Please read carefully the description of each configuration variable 
and modify them as needed. 


### Run RACA 

There is a wrapper Perl script, `Run_RACA.pl`. To run RACA, type as:

    <path to RACA>/Run_RACA.pl <path to the configuration file> 
 
Example dataset (Tibetan antelope assembly)
-------------------------------------------

### Download the dataset

Visit http://bioinfo.konkuk.ac.kr/RACA/ and click the 
link "Tibetan antelope (TA) data". Then you can download the file 
`TAdata.tgz` file. 

### Compile the dataset

Go into the directory where you downloaded the TAdata.tgz file and run:

    tar xvfz TAdata.tgz
    cd TAdata/
    make

### Run RACA for the dataset

In the TAdata directory run:

    <path to RACA>/Run_RACA.pl params.txt

### Where are output files?

In the `TAdata/Out_RACA` directory.
		

What are produced?
------------------

In the output directory that is specified in the above configuration file, 
the following files are produced.

- `rec_chrs.refined.txt `

    This file contains the order and orientation of target scaffolds in 
    each reconstructed RACA chromosome fragment. Each column is defined 
    as:  

    +  Column1: the RACA chromosome fragment id
    +  Column2: start position (0-based) in the RACA chromosome fragment
    +  Column3: end position (1-based) in the RACA chromosome fragment 
    +  Column4: target scaffold id or 'GAPS'
    +  Column5: start position (0-based) in the target scaffold
    +  Column6: end position (1-based) in the target scaffold

- `rec_chrs.<ref_spc>.segments.refined.txt`

    This file contains the mapping between the RACA chromosome fragments 
    and the genome sequences of the reference species `<ref_spc>`. 

- `ref_chrs.<tar_spc>.segments.refined.txt`
    
    This file contains the mapping between the RACA chromosome fragments 
    and the genome sequences of the target species `<tar_spc>`. 

- `ref_chrs.<out_spc>.segments.refined.txt`
    
    This file contains the mapping between the RACA chromosome fragments 
    and the genome sequences of the outgroup species `<out_spc>`. This file 
    is created for each outgroup species. 

- `rec_chrs.adjscores.txt`

    This file constins the adjacency scores that were used to reconstruct 
    the RACA chromosome fragments. Each column is defined as:

    +  Column1: the RACA chromosome fragment id
    +  Column2: start position (1-based) in the RACA chromosome fragment
    +  Column3: end position (1-based) in the RACA chromosome fragment 
    +  Column4: the adjacency score

- `rec_chrs.size.txt`

    This file contains the total size (the second column) and the total 
    number of target scaffolds (the third column) that are placed in each 
    RACA chromosome fragment (the first column).  

There are other intermediate files and directories in the output directory.
They can be safely ignored.  

How to ask questions?
---------------------

Contact Jaebum Kim (Jaebum.Kim@gmail.com)

Citation
--------

Kim J, Larkin DM, Cai Q, Asan, Zhang Y, Ge RL, Auvil L, Capitanu B, Zhang G,
Lewin HA, Ma J. (2013) [Reference-assisted chromosome assembly][PNAS].
 _Proc Natl Acad Sci USA_ **110**(5):1785-90.<br/>
doi: 10.1073/pnas.1220349110; PMID: 23307812; PMCID: PMC3562798.

[PNAS]: https://www.pnas.org/content/110/5/1785
