use warnings;
use diagnostics;
if(@ARGV != 3){
	die "*.pl <kssdcomposite_output.tsv> <best.gtdbr207_psid2ncbi_specid.tsv> <ncbitaxid_rank_parentnode_name.gtdbr207_pseudoidrelated.tsv>";
}
#Ranks:superkingdom|phylum|class|order|family|genus|species
$median_thr = 1;
$avgpct9899_thr = 3;
$shkm_thr = 7;
#LOW THRESHOULD
$low_avgpct9899_thr = 2;
$small_val = 0.001;

open $nodef, $ARGV[2] || die "cant open $ARGV[2]: $!";
while(<$nodef>){
	chomp;
	($node,$rank,$pa,$name)=(split /\t+/)[0..3];
	$node2rank{$node} = $rank;
	$node2pa{$node} = $pa;
	$node2name{$node} = $name;
#	push @{$rank_cate{$rank}},$node;
}
close $nodef;

open $spcidf, $ARGV[1] || die "cant open $ARGV[1]: $!";
while(<$spcidf>){
	chomp;
	($psid,$ncbi_id)=(split /\t+/)[0,1];
	$psid2ncbi_sp{$psid} = $ncbi_id;
}
close $spcidf;

#use origin kssd composite output
#format:Qry\tRef\tShare_kmc\tAvg_kmc\tAvgpct9899_kmc\tMedian_kmc\tPct98_kmc
open $kssdf, $ARGV[0] || die "cant open $ARGV[0]: $!";

while(<$kssdf>){
	chomp;
	($sample,$ref, $shkm, $avgpct9899,$median)= (split /\t+/)[0,1,2,4,5];
	$sample=~s/[^0-9a-zA-Z_.]/_/g;
	$psid=(split /_/, $ref)[0];
	if($avgpct9899 > $avgpct9899_thr and $median > $median_thr and $shkm > $shkm_thr){
		$data{$sample}->{$psid} = $avgpct9899 - $avgpct9899_thr;
		$sum{$sample} += ($avgpct9899 - $avgpct9899_thr);
	}
#assign small value to low abundance species to increase completeness
	elsif($avgpct9899 >= $low_avgpct9899_thr and $shkm > $shkm_thr ){
		$depth = $avgpct9899 - $avgpct9899_thr > $small_val ? $avgpct9899 - $avgpct9899_thr : $small_val;
		$data{$sample}->{$psid} = $depth;		
		$sum{$sample} += $depth; 
	}
}

close $kssdf;

@ranks = ("superkingdom","phylum","class","order","family","genus","species");
foreach $rank(@ranks){$rankidx{$rank} = 1;};

foreach $sample(keys %data){
	foreach $rank (@ranks){
		 $rank_cate{$rank} = ();
	} ;
	%nctax_ab = ();
	foreach $psid( keys %{$data{$sample}}){
		$nc_spcid = $psid2ncbi_sp{$psid};
		push @{$rank_cate{$node2rank{$nc_spcid}}},$nc_spcid if !exists $nctax_ab{$nc_spcid};
		$nctax_ab{$nc_spcid} = $data{$sample}->{$psid} / $sum{$sample} * 100 ;
		
		$node = $node2pa{$nc_spcid};
		while($node != 1){			
			push @{$rank_cate{$node2rank{$node}}},$node if !exists $nctax_ab{$node};
			$nctax_ab{$node} += $nctax_ab{$nc_spcid};
			$node = $node2pa{$node};	
		}
	}

#profiling output

print "# Taxonomic Profiling Output\n";
print "\@SampleID:",$sample,"\n";
print "\@Version:0.9.1\n";
print "\@Ranks:superkingdom|phylum|class|order|family|genus|species\n";
print "\@TaxonomyID:ncbi-taxonomy_2021.07.19\n";
print "\@__program__:kssd2\n";
print "\@\@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n";		
	foreach $ele (@ranks){
		  @sorted_taxid = sort{$nctax_ab{$b} <=> $nctax_ab{$a} } @{$rank_cate{$ele}} ;
			foreach $taxid (@sorted_taxid) {
				$node = $taxid;
				@tmp_arr = ();
				@tmp_arr_name = ();
				while($node2pa{$node} != 1 ){
					if ( exists $rankidx{$node2rank{$node}} ) {
						push @tmp_arr,$node;
						push @tmp_arr_name, $node2name{$node};
					} 
					$node = $node2pa{$node};
				};
				$taxpath = join '|', reverse @tmp_arr;
				$taxpathsn = join '|', reverse @tmp_arr_name;
				print  $taxid,"\t",$ele,"\t",$taxpath,"\t",$taxpathsn,"\t",sprintf("%.4f", $nctax_ab{$taxid}),"\n";
			
			}
	}
}































