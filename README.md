# Kssd: K-mer substring space sampling/shuffling Decomposition

Kssd is a command-line tool for large-scale sequences sketching and resemblance- and containment-analysis. It sketches sequences by k-mer substring space sampling/shuffling. It handles DNA sequences of both fasta or fastq format, whether gzipped or not. Kssd run on linux system, currently not support for MacOS and Windows OS.
### Table of contents
1.  [Installation](#1-installation)
2.  [Quick Tutorial](#2-quick-tutorial)
3.  [For Advanced Users](#3-for-advanced-users)
    1. [K-mer substring space shuffling](#31-k-mer-substring-space-shuffling)
    2. [Sketching sequences](#32-sketching-sequences)
        1.  [Sketching references](#321-sketching-references)  
        2.  [Sketching queries](#322-sketching-queries)
    3.  [Distance estimation](#33-distance-estimation)
        1.  [Reference against references distance](#331-reference-against-references-distance) 
        2.  [Search the queries against the references](#332-search-the-queries-against-the-references)
    4.	[Combine Queries](#34-combine-queries)
4.  [How to cite](#4-how-to-cite)     

# 1 Installation 
```
git clone git@github.com:yhg926/public_kssd.git &&
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

Yi, H., Y. Lin, and W. Jin, Sequences Dimensionality-Reduction by K-mer Substring Space Sampling Enables Effective Resemblance- and Containment-Analysis for Large-Scale omics-data. bioRxiv, 2019: p. 729665. (https://www.biorxiv.org/content/10.1101/729665v1)


