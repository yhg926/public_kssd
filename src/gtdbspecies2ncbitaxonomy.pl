use warnings;
use diagnostics;
if(@ARGV != 2){
	die " *.pl <gtdb_species_list> <gtdb2ncbi_taxonomy.tsv>";
}

open $db,$ARGV[1] || die "cant open $ARGV[1] : $!";
while(<$db>){
	chomp;
	($gtdb,$ncbi_id,$ncbi_tax)=(split /\t+/)[0,1,2];

	if($gtdb =~ /\;s__(.+)/ ){
		if (exists $hash{$1}) {
			$hash{$1} = $hash{$1}."|".$ncbi_id."_".$ncbi_tax ;
		}
		else{
			$hash{$1} = $ncbi_id."_".$ncbi_tax ;
		}
	}
}
close $db;

open $f,$ARGV[0] || die "cant open $ARGV[0] : $!";

while($gtdb_sp=<$f>){
	chomp $gtdb_sp;
	print $gtdb_sp,"\t";
	if(exists $hash{$gtdb_sp}) {
		print $hash{$gtdb_sp};
	}
	else {
		print "0";
	}
	print "\n";
}
close $f;
