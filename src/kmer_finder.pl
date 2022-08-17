use warnings;
use diagnostics;
if(@ARGV != 2){
	die "*.pl <*.fasta> <K> ";
}
open $fas, $ARGV[0] || die "can't open $ARGV[0] :$!";
$k=$ARGV[1];

$/ = '>';
while( $line=<$fas>){
	chomp $line;
	@tmp = split('\n',$line) ; 
	shift @tmp;
	$read = join '',@tmp;
	for($i=0;$i< length($read) - $k+1 ;$i++){
		$kmer = substr($read,$i,$k);
		$rckmer = $kmer;
		$rckmer  =~ tr/ACGTacgt/TGCAtgca/;		
		$rckmer = reverse $rckmer ;
		$unikmer = $kmer lt $rckmer?$kmer:$rckmer;
		$hash{$unikmer} = 1; 
	}		
}
foreach $ele (keys %hash) {
	print $ele,"\n";
}	
