use warnings;
use diagnostics;
if(@ARGV!=2){
	die "*.pl <genomes.gtdbpseudotaxid2ncbitax.tsv> < first(0) or all(1)>";
}
open $gs,$ARGV[0] || die "cant open $ARGV[0]:$!";
while(<$gs>){
	chomp;
	next if /GTDB_AC/;

	($gid,$psid,$gtname,$ncbi_tids,$ncbi_taxn)=(split /\t+/)[0,1,2,3,4];	
  next if $gid !~ /GC[AF]_\d+/ ;
	$gtdbkey=$psid."_".$gtname;

	@tids = (split /\|/,$ncbi_tids);
	$len = @tids - 1 ;
	$ncbi_tid_path = join '|',@tids[0..$len-1];
	$ncbi_spcid = $tids[$len-1];
	
	@tnames = (split /\|/,$ncbi_taxn);
#	$len = @tnames <= 6 ? @tnames - 1 : 6; 
	$ncbi_tname_path = join '|',@tnames[0..$len-1]; 
	
	$hash{$gtdbkey}->{$ncbi_spcid}->{"gn"}++;
	$hash{$gtdbkey}->{$ncbi_spcid}->{"tid_path"} = $ncbi_tid_path ;
	$hash{$gtdbkey}->{$ncbi_spcid}->{"tname_path"} = $ncbi_tname_path; 	
}
close $gs;

foreach $gtdbkey (sort keys %hash){
	@sorted_spcid	= sort{ $hash{$gtdbkey}->{$b}->{"gn"} 
			 <=> $hash{$gtdbkey}->{$a}->{"gn"} } keys %{$hash{$gtdbkey}};

	if ($ARGV[1] >0){
		foreach $ele (@sorted_spcid){
			print $gtdbkey,"\t",$hash{$gtdbkey}->{$ele}->{"tid_path"},"\t",
				 $hash{$gtdbkey}->{$ele}->{"tname_path"},"\t", $hash{$gtdbkey}->{$ele}->{"gn"},"\n";
		}
	}
	else{
		 print $gtdbkey,"\t",$hash{$gtdbkey}->{$sorted_spcid[0]}->{"tid_path"},"\t",
         $hash{$gtdbkey}->{$sorted_spcid[0]}->{"tname_path"},"\t", $hash{$gtdbkey}->{$sorted_spcid[0]}->{"gn"},"\n";
	}
	
}
