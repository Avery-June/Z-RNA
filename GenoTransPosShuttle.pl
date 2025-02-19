#!/usr/bin/perl
use 5.010;
use warnings;
use strict;
use diagnostics;
use Getopt::Long;

my ($opt,$gpe,$bin,$in);
my $typeofID = 1;
GetOptions(
    "opt|o:s"       =>      \$opt,
    "gpe|g:s"       =>      \$gpe,
    "bin|b"         =>      \$bin,
    "in|i:s"        =>      \$in,
    "idtype|t:i"    =>      \$typeofID,
    "help|h"        =>      sub{&usage;exit(-1);}
);

if (!defined $opt || ($opt ne "trans2geno" && $opt ne "geno2trans")){
    print STDERR "please set --opt|-o by 'trans2geno' or 'geno2trans'.\n";
    exit(-1);
}
if($typeofID != 1 && $typeofID != 2){
    print STDERR "please set --idtype|t, 1 for trans, 2 for geno|trans.\n";
    exit(-1);
}
if(!defined $gpe){exit(-1);}

open GPE,$gpe;

my %struc;
while(<GPE>){
	chomp;
	next if (/^#/);
	my @struc = split("\t",$_);
        if($bin){@struc = @struc[1..$#struc];}
	my ($enst,$chrom,$strand,$txStart,$txEnd,$exonCount,$exonStarts,$exonEnds,$ensg) = @struc[0..4,7..9,11];

	my @starts = split(",",$exonStarts);
	my @ends = split(",",$exonEnds);

        my $id;
        if ($typeofID == 2){$id = "$ensg|$enst";}
        if ($typeofID == 1){$id = "$enst";}

	$struc{"$id"}->{'chrom'} = $chrom;
	$struc{"$id"}->{'strand'} = $strand;
	$struc{"$id"}->{'txStart'} = $txStart;
	$struc{"$id"}->{'txEnd'} = $txEnd;
	$struc{"$id"}->{'exonCount'} = $exonCount;

	$struc{"$id"}->{'starts'} = \@starts;
	$struc{"$id"}->{'ends'} = \@ends;
# No. of exon (1-base)
	for(my $i = 1;$i <= $exonCount;$i++){
	    $struc{"$id"}->{'exon'}->{$i}->{'start'} = $starts[$i-1];
	    $struc{"$id"}->{'exon'}->{$i}->{'end'} = $ends[$i-1];
	}

}

my $IN;
unless(defined $in){$IN = \*STDIN;}
else{open $IN,$in;}

if($opt eq 'trans2geno'){
    while(<$IN>){
        chomp;
        next if (/^#/);
        my ($id,$pos) = split;
        if($struc{$id}->{'strand'} eq "-"){
            $pos = &length($struc{"$id"}->{'starts'},$struc{"$id"}->{'ends'}) - $pos + 1;
        }
        print "$_\t";
        print $struc{"$id"}->{'chrom'};
        print "\t";
        print &trans2geno($id,$pos);
        print "\n";
    }
}

if($opt eq 'geno2trans'){
    while(<$IN>){
        chomp;
        next if (/^#/);
        my ($id,$chr,$pos) = split;
        print "$_\t";
        print &geno2trans($id,$pos);
        print "\n";
    }
}

# FUNCTION: \@start = (1,2,3) \@end = (2,3,4)
#           return each(@end-@start)
sub length{
	my ($s,$e)=(@_);
	my $total_length;
	for (my $k = 0;;$k++){
		if (defined($s->[$k])){$total_length += $e->[$k]-$s->[$k];}
		else {last;}
	}
	return $total_length;
}

# FUNCTION: IN: transID, transPos; OUT: genoPos
sub trans2geno{
	my ($RNAME,$POS) = @_;
	my ($LENG,$SUM) = (0,0);
	for (my $j=1;;$j++){
		$LENG = $struc{"$RNAME"}->{'exon'}->{$j}->{'end'} - $struc{"$RNAME"}->{'exon'}->{$j}->{'start'};
		$SUM += $LENG;
		if ($POS <= $SUM){
			return ($struc{"$RNAME"}->{'exon'}->{$j}->{'start'} + $POS - ($SUM - $LENG));
                        last;
		}
	}
}

# FUNCTION: IN: transID,genoPos; OUT: transID
sub geno2trans{
    my ($RNAME,$POS) = @_;
    my ($LENG,$SUM) = (0,0);
    for(my $i = 1;;$i++){
        $LENG = $struc{$RNAME}->{'exon'}->{$i}->{'end'} - $struc{$RNAME}->{'exon'}->{$i}->{'start'};
        $SUM += $LENG;
        if ($POS > $struc{$RNAME}->{'exon'}->{$i}->{'start'} && $POS <= $struc{$RNAME}->{'exon'}->{$i}->{'end'}){
            if($struc{$RNAME}->{'strand'} eq "-"){
                return(&length($struc{$RNAME}->{'starts'},$struc{$RNAME}->{'ends'})-($SUM -$LENG + $POS - $struc{$RNAME}->{'exon'}->{$i}->{'start'})+1);
            }else{
                return($SUM - $LENG + $POS - $struc{$RNAME}->{'exon'}->{$i}->{'start'});
            }
            last;
        }
    }
}


# FUNCTION: Usage
sub usage{
        print STDERR <<HD_USAGE
[Usage]:
        $0 --opt [] --gpe .gpe (-b) --in .pos -t [1] --help
[parameters]:
    --opt|-o:s              trans2geno or geno2trans?
    --gpe|-g:s              .gpe gene structure
    --bin|-b                with bin column?
    --in|-i:s               transpos (transID,1-base position); genopos (transID,chr,1-base postion)
    --idtype|-t:i           type of ID (1 for 'trans', 2 for 'gene|trans') [1]
    --help|-h               print this help message
[Release Notes]:
	\@version 1.0 (01-05-2015)
		beta version
HD_USAGE
}
