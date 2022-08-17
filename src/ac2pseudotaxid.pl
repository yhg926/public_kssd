use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <accessions> <ac2pseudotaxid.gtdb310k_bac_ar>";
}

open $tidf,$ARGV[1]||die "cant open $ARGV[1]:$!";
while(<$tidf>){
	chomp;	
	($ac,$tid,$name)=(split /\t+/)[0,1,2];
	if(!defined $name){
		$hash{$ac} = $tid;
	}
	else {
		$hash{$ac} = $tid."\t".$name;
	}
}
close($tidf);

open $acf,$ARGV[0]||die "cant open $ARGV[0]:$!";

while($ac=<$acf>){
	chomp $ac;
	print $ac,"\t";
	if(exists $hash{$ac}){
		print $hash{$ac};
	}
	else{
		print "0";
	}
	print "\n";

}
close($acf);
