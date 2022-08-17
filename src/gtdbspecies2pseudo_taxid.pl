use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <gtdb species name> <pseudotaxid2species.tsv>";
}

open $tidf,$ARGV[1]||die "cant open $ARGV[1]:$!";
while(<$tidf>){
	chomp;
	($tid,$name)=(split /\t+/)[0,1];
	$tid=~s/\s+//g;
	$hash{$name} = $tid;
}
close($tidf);

open $nf,$ARGV[0]||die "cant open $ARGV[0]:$!";

while(<$nf>){
	chomp;
	$name=(split /\t+/)[0];
	if(exists $hash{$name}){
		print $hash{$name};
	}
	else{print "0";}
	print "\t$name\n";
}
close($nf);
