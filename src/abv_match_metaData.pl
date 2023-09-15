use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <abv search output> <All_metag.run_proj_sample_orgnism.cs>";
}

open $meta,$ARGV[1] || die "cant open $ARGV[1]:$!";
$head=<$meta>;
chomp $head;
($colname2,$colname3,$colname4)=(split /\,/,$head)[1,2,3];
while(<$meta>){
	chomp;
	($run,$bpj,$bsmp,$orgnism) = (split /\,/);
	if((defined $bsmp) and (defined $orgnism)) {
		$hash{$run} = $bpj."\t".$bsmp."\t".$orgnism;
	}
}
close($meta);

open $abv,$ARGV[0] || die "cant open $ARGV[0]:$!";

while(<$abv>){
	chomp;
	($abvname,$measure)=(split /\t/);
	if ($abvname !~ /\.abv$/){
		print $abvname,"\t",$measure,"\t",$colname2,"\t",$colname3,"\t",$colname4,"\n";
	}
	else{
		$abvname=~s/\.abv$//;
		print $abvname,"\t",$measure,"\t";
		if(!exists $hash{$abvname}){
			print "NA\tNA\tNA\n";
		}
		else {
			print $hash{$abvname},"\n";
		}
	
	}	
}
close($abv);



