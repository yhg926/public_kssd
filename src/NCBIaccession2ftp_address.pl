use warnings;
use diagnostics;
if(@ARGV!=1){
	die "*.pl <accession_ASMid.tsv>";
}

open $f,$ARGV[0] || die "cant open $ARGV[0] :$!";

while(<$f>){
	chomp;
	($ac,$asm )=(split /\t+/)[0,1];
	($fac)=($ac=~/(GC[AF]_[0-9.]+)/);
	$asm=~s/\s/_/g;
	($gc,$num)=(split /_/,$fac)[0,1];
	($n1,$n2,$n3)=($num=~/(\d{3})(\d{3})(\d{3})/);
		

	print "rsync://ftp.ncbi.nlm.nih.gov/genomes/all/",$gc,'/',$n1,'/',$n2,'/',$n3,'/',$fac,'_',$asm,'/',$fac,'_',$asm,"_genomic.fna.gz","\n";	

}
