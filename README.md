# Kssd: K-mer substring space sampling/shuffling Decomposition

Kssd is a command-line tool for large-scale sequences sketching and resemblance- and containment-analysis. It sketches sequences by k-mer substring space sampling/shuffling, please see Methods part of our Genome Biology paper (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02303-4) for how it works. It handles DNA sequences of both fasta or fastq format, whether gzipped or not. Kssd run on linux system, currently not support for MacOS and Windows OS.
### Table of contents
1.  [Installation](#1-installation)
2.  [Quick Tutorial](#2-quick-tutorial)
3.  [For Advanced Users](#3-for-advanced-users)
    1. [K-mer substring space shuffling](#31-k-mer-substring-space-shuffling)
    2. [Sketching sequences](#32-sketching-sequences)
        1.  [Sketching references](#321-sketching-references)  
        2.  [Sketching queries](#322-sketching-queries)
	3.  [Sketching from data streaming](#323-Sketching-from-data-streaming)
	4.  [Set operations](#324-Set-operations)
		1.  [Sketches union](#3241-Sketches-union)
		2.  [Sketches intersection](#3242-Sketches-intersection)
		3.  [Sketches subtraction](#3243-Sketches-subtraction)	 
    3.  [Distance estimation](#33-distance-estimation)
        1.  [Reference against references distance](#331-reference-against-references-distance) 
        2.  [Search the queries against the references](#332-search-the-queries-against-the-references)
    4.	[Combine Queries](#34-combine-queries)
4.  [How to cite](#4-how-to-cite)     

# 1 Installation 
```
git clone https://github.com/yhg926/public_kssd.git &&
cd public_kssd &&
make 
```
# 2 Quick-Tutorial
```
cd test_fna;
#sketch references
 ../kssd dist -L 3 -r ./seqs1 -o reference
#sketch queries using same kmer substring shuffled space with the references (.shuf file)
../kssd dist -L reference/default.shuf -o query ./seqs2
#Search queries against references db 
../kssd dist -r reference/ref -o distout query/qry
# or you can compute the pairwise distance of references
../kssd dist -r reference/ref -o reference reference/qry
```
Here is the explanation of the output file "distance.out" (please see [How to cite](#4-how-to-cite) for the referred equations)

Column | Explanation
---|---
Qry | query
Ref  | reference
Shared_k\|Ref_s\|Qry_s | number of shared k-mer between the sketches of the reference and the query\|reference sketch-size\|query sketch-size 
Jaccard|Jaccard-coefficient (Eq. 2)
MashD| mash distance  (Eq. 4)
ContainmentM| containment-measurement(Eq. 3)
AafD| aaf-distance (Eq. 5) 
Jaccard_CI| 0.95 confidence intervel for Jaccard-coefficient
MashD_CI| 0.95 confidence intervel for mash distance
ContainmentM_CI| 0.95 confidence intervel for containment-measurement
AafD_CI| 0.95 confidence intervel for aaf-distance
P-value(J)| p-value for Jaccard-coefficient(Eq. 14)
P-value(C)| p-value for containment-measurement(Eq. 14)
FDR(J)| false discover rate for Jaccard-coefficient
FDR(C)| false discover rate for containment-measurement

# 3 For Advanced Users
## 3.1 K-mer substring space shuffling
```
kssd shuffle -k <half_length_of_k-mer> -s <half_length_of_k-mer_substring> -l <dimensionality-reduction_level > -o <shuffled_k-mer_substring_space_file>
```
This step can be omitted, and you can skip to step 2 if you wish to use default setting of `-s`. Other wise read below:  
This command will generate a file suffixed by ‘.shuf’ which keeps the shuffled k-mer substring space, this file would then took as input for sequences sketching or decomposition.  
`-k`: Half-length of k-mer, `-k x` meaning use k-mer of length `2x`. For bacterial `-k 8` is recommand; for mammals, `-k 10` is recommand; for other genome size in-between, `-k 9` is recommand.  
`-s`: Half-length of k-mer substring, `-s x` meaning the whole space is the collection of all `2x-mer`. Make sure `l < s < k`. The default setting is `-s 6`, usually there is no need to change this setting.   
`-l`: The level of dimensionality-reduction. `-l x` meaning the expected rate of dimensionality-reduction is `$16^x$`; for bacterial `-l 3` is recommand; for mammals, `-l 4` or `-l 5` is recommand. `l < s`.  
`-o` output .shuf file.
## 3.2 Sketching sequences
### 3.2.1 Sketching references
```
kssd dist -r <.fasta/fastq_dir> -L <.shuf_file or dimentionality-reduction_level> [-k <half_k-mer_length>] -o <outdir>
```
`-L`: The`.shuf` file generated from `kssd shuffle` or the the level of dimensionality-reduction.  
 
  If you feed `-L` the `.shuf` file, there is no need to specify `-k` again, since it has already been set in the `.shuf` file.
  
  Else if you feed `-L` the level of dimensionality-reduction, new `.shuf` file will generated and used. Actually, command:
```
kssd dist -r <.fasta/fastq_dir> -L <dimentionality-reduction_level> -k <half_length_of_k-mer> -o <ref_outdir>
```
is equivalent to  
```
kssd shuffle -k <half_length_of_k-mer> -s <half_length_of_k-mer_substring> -l <dimensionality-reduction_level> -o <ref_outdir/default.shuf> &&
kssd dist -r <.fasta/fastq_dir> -L <ref_outdir/default.shuf> -o <ref_outdir>
```
The expected rate of dimensionality-reduction for `-L x` is `$16^x$`; for bacterial `-L 3` is recommand; for mammals, `-L 4` or `-L 5` is recommand. 

`-r`: Feed it with the sequences (fasta or fastq, gzipped or not) that you want built as the references-db.  
  
`-k`: Half-length of k-mer, `-k x` meaning use k-mer of length `2x`. For bacterial `-k 8` is recommand; for mammals, `-k 10` is recommand; for other genome size in-between, `-k 9` is recommand.  
  
`-o`: There are two folders `ref/` and `qry/` in the output dir `ref_outdir`.  In Step 3 distance estimation `ref_outdir/ref` feed as references for `-r` and `ref_outdir/qry` feed as queries   

### 3.2.2 Sketching queries
To compare queries with references, queries need be skeched using the same `.shuf` file with that of references.
```
kssd dist -o <qry_outdir> -L <ref_outdir/default.shuf or the_.shuf_file_used_by_references> <queries_.fasta/fastq_dir>
```
`-o`: There is only one folder `qry/` in the output dir `qry_outdir`. In Step 3 distance estimation `qry_outdir/qry` feed as queries.

### 3.2.3 Sketching from data streaming
Suppose you want sketching Sequence Read Archive Accesssion ERR000001, just run:
```
kssd dist -L <your .shuf file> -n 2 -o <outdir> --pipecmd "fastq-dump --skip-technical --split-spot -Z" ERR000001
``` 
or prefetch first
```
prefetch ERR000001 && kssd dist -L <your .shuf file> -n 2 -o <outdir> --pipecmd "fastq-dump --skip-technical --split-spot -Z" <.sra dir>/ERR000001.sra
```
### 3.2.4 Set operations
#### 3.2.4.1 Sketches union 
```
kssd set -u <qry_outdir/qry> -o <union_outdir>  
```
It will create the union sketch in <union_outdir> from the combined queries sketch in <qry_outdir/qry>. Note the combined queries sketch is just a sketch combined from all queries sketches, the union operation deduplicate those integers duplicated in different queries;
#### 3.2.4.2 Sketches intersection
```
kssd set -i <union_outdir> -o <intersect_outdir> <qry_outdir/qry>
```
It will create the intersection sketch in <intersect_outdir> between the union sketch in <union_outdir> and the combined queries sketch in <qry_outdir/qry>;
#### 3.2.4.3 Sketches subtraction
```
kssd set -s <union_outdir> -o <subtract_outdir> <qry_outdir/qry>
```
It subtracts the union sketch in <union_outdir> from the combined queries sketch in <qry_outdir/qry> and creates the remainder sketch in <subtract_outdir>

## 3.3 Distance estimation
### 3.3.1 Reference against references distance
If you only want to compute pairwise distances of all references, run:
```
kssd dist -r <ref_outdir/ref> -o <outdir> <ref_outdir/qry>
```
### 3.3.2 Search the queries against the references
Or if you want search the queries against the references, run:
```
kssd dist -r <ref_outdir/ref> -o <outdir> <qry_outdir/qry>
```
The `ref_outdir` and `qry_outdir` are the sketches created in step 2.  
   
  
  The distance will output to `<outdir>/disntance`

## 3.4 Combine Queries
If you have queries generated from different running batches, you can combine them by:
```
kssd dist -o <outdir> <path_to_query_batch1> [<path_to_query_batch2> ...]  
```
Make sure all queries batches use the same .shuf file

# 4. How to cite

Yi, H., Lin, Y., Lin, C. et al. Kssd: sequence dimensionality reduction by k-mer substring space sampling enables real-time large-scale datasets analysis. Genome Biol 22, 84 (2021). https://doi.org/10.1186/s13059-021-02303-4



