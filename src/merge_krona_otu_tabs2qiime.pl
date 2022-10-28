use warnings;
use diagnostics;
if(@ARGV < 2 || (-f $ARGV[0])){
	die "*.pl <outdir> <Krona_in.tsv> [<Krona_in.tsv> ...]";
}
elsif(! -d $ARGV[0]){
	mkdir($ARGV[0]);
}

$otuid = 0;

for ($i=1; $i< @ARGV; $i++) {
	open $fp, $ARGV[$i] || die "cannot open $ARGV[$i]:$!";

		while(<$fp>){
			chomp;
			@tmp_row = (split /\t/);
			$val  = shift @tmp_row ;
			$taxa = join ';',@tmp_row ;
	
			if (!exists $taxa2otuid{$taxa}){
				$taxa2otuid{$taxa} = $otuid;
				push @otuid2taxa,$taxa;
				$otuid++;
			}
			$abund[$taxa2otuid{$taxa}][$i-1] = $val ;
			
		}
	close $fp;
}

open OTU, ">$ARGV[0]/otu.tsv"; 
open TAX, ">$ARGV[0]/taxonomy.tsv";
open META, ">$ARGV[0]/meta.tsv";
print META "sample-id","\n";

print OTU "#OTU";
for ($i=1;$i<@ARGV; $i++){
	$tempname= $ARGV[$i] ;
	$tempname =~ s/.*\///g;
	$tempname =~ s/\..*//g;
	print OTU "\t", $ARGV[$i];
	print META $ARGV[$i],"\n";
}
print OTU "\n";

for ($i=0;$i<$otuid; $i++){
	print OTU "OTU_",$i;
	print TAX "OTU_",$i,"\t",$otuid2taxa[$i],"\n";

	for ($j=1;$j<@ARGV; $j++){
		print OTU "\t",defined $abund[$i][$j-1] ? $abund[$i][$j-1] : 0; 
	}
	print OTU "\n";
}
close OTU;

close TAX;
close META;














