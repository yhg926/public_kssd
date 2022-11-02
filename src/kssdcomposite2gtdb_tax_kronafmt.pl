use warnings;
use diagnostics;
if(@ARGV != 2){
	die "*.pl <kssdcomposite_output.tsv> <gtdbr207_psid2krona_taxonomy.tsv>";
}
#Ranks:superkingdom|phylum|class|order|family|genus|species
$median_thr = 1;
$avgpct9899_thr = 3;
$shkm_thr = 8;
#LOW THRESHOULD
$low_avgpct9899_thr = 2;
$small_val = 0.1;

open $taxf,$ARGV[1] || die "can't open $ARGV[1] :$!";

while(<$taxf>){
	chomp;
	@line =(split /\t+/);
	$psid = shift @line ;
	$psid=~s/ //g; 
	
	$hash{$psid} = join "\t",@line;

}

close $taxf;

open $kssdf,$ARGV[0] || die "can't open $ARGV[0] :$!";

use File::Basename;
@cmp_fmt = (".gz"); #or other compress format
@seq_fmt = (".fq",".fastq",".fa",".fna",".fas",".fasta");
$num_smp=0;
$sample = "NULL";
while(<$kssdf>){
        chomp;
        ($sample,$ref, $shkm, $avgpct9899,$median)= (split /\t+/)[0,1,2,4,5];
				$sample = basename($sample,@cmp_fmt);
			  $sample = basename($sample,@seq_fmt);

				if (!exists $dic{$sample}) {
					$dic{$sample} = 1; 
					$num_smp++;
					if ($num_smp >1){
          	printf("Error: Client mode only accept 1 sample one time\n");
          	exit(1);
        	}					
				}

        $psid=(split /_/, $ref)[0];
				next if $shkm <= $shkm_thr; 
        if($avgpct9899 > $avgpct9899_thr and $median > $median_thr){
                $depth{$psid} = $avgpct9899 - $avgpct9899_thr;
                $sum += $depth{$psid};
        }
#assign small value to low abundance species to increase completeness
        elsif($avgpct9899 >= $low_avgpct9899_thr){
                $depth{$psid} = $avgpct9899 - $avgpct9899_thr > $small_val ? $avgpct9899 - $avgpct9899_thr : $small_val;               
                $sum += $depth{$psid};
        }
}

close $kssdf;

$dirname = "./temp_Krona_fmt";
mkdir $dirname, 0755;
open $outf, "> $dirname/$sample";

foreach $psid ( sort{ $depth{$b} <=> $depth{$a} or $a <=> $b  } keys %depth ){	
	printf $outf "%.4f\t%s\n",  $depth{$psid}/$sum , $hash{$psid};  
}

close $outf;
