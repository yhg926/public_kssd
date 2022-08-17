#!/usr/bin/perl
use warnings;
use diagnostics;
if(@ARGV!=2) {
	die "*.pl <genomelist> <all_genome2taxid_file.tsv>";
}

open $db,$ARGV[1]||die "can't open $ARGV[1]:$!";
while(<$db>){
	chomp;
	($gid,$taxid,$taxname)=(split /\t+/)[0,1,2];

	if (defined $taxname){
		$str = $taxid."\t".$taxname;
	}
	else{
		$str = $taxid;
	}		
	$hash{$gid} = $str;
		
}
close($db);

open $g,$ARGV[0]||die "can't open $ARGV[0]:$!";

while($rawid=<$g>){
	chomp $rawid;
	($gid) = ($rawid =~ /(GC[AF]_[0-9.]+)/) ;
	$printstr = exists $hash{$gid} ? $hash{$gid} : 0;
	print $gid,"\t",$printstr,"\n";
}
close($g);
