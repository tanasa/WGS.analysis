#!/usr/bin/perl
#use strict;

# filtering function of AD and AF, SOR and FS.

use List::Util qw[min max];
use warnings;
use strict;

###########################################
###########################################

my $inputfile = '';
my $outputfile = '';
my $outputfile2 = '';

my $line = '';
my @lines = ();

my @CHROM = ();
my @POS = ();
my @ID = ();
my @REF = ();
my @ALT = ();
my @QUAL = ();
my @FILTER = ();
my @INFO = ();
my @FORMAT = ();
my @TUMOR = ();
my @NORMAL = ();

my $MIN_RD = 0;
my $MAX_RD = 0;

my $MAX_FS = 0;
my $MAX_SOR = 0;

my $MIN_AD = 0;
my $MIN_AF = 0;

#####################################
#####################################

my $number_SNP_unfiltered = 0;
my $number_INSERTION_unfiltered = 0;
my $number_DELETION_unfiltered = 0;
my $number_VARIANTS_unfiltered = 0;

my $number_SNP_filtered = 0;
my $number_INSERTION_filtered = 0;
my $number_DELETION_filtered = 0;
my $number_VARIANTS_filtered = 0;

#####################################
### reading the INPUT file (VCF file)

$inputfile = $ARGV[0];

#$MIN_RD = $ARGV[1];
#$MAX_RD = $ARGV[2];
#$MAX_FS = $ARGV[3];
#$MAX_SOR = $ARGV[4];
#$MIN_AD = $ARGV[5];
#$MIN_AF = $ARGV[6];

$MIN_RD = 10;
$MAX_RD = 1000;

$MAX_FS = 60;
$MAX_SOR = 4;

$MIN_AD = 5;
$MIN_AF = 0.05;

##################################
##################################

#$outputfile1 = $inputfile.".filtered.FS.SOR.vcf";
#$outputfile2 = $inputfile.".filtered.FS.SOR.RD.AD.AF.vcf";

$outputfile = $inputfile.".perl.filtering.FS.SOR.and.AD.AF.vcf";
$outputfile2 = $inputfile.".perl.filtering.FS.SOR.and.AD.AF.summary";

open(INPUT, "< $inputfile") or die "cannot open $inputfile : $!\n";
open(OUTPUT, "> $outputfile");
open(OUTPUT2, "> $outputfile2");

#open(OUTPUT1, "> $outputfile1");
#open(OUTPUT2, "> $outputfile2");

###############################################
###############################################
### reading the parameters to filter:
### MIN_RD
### MAX_RD
### FS
### SOR
### MIN_AD
### MIN_AF
###############################################
###############################################

my @t = ();
my @v = ();
my @INFO_array_length = ();
my @INFO_array = ();

my @FORMAT_array_length = ();
my @FORMAT_array = ();

my @NORMAL_array_length = ();
my @NORMAL_array = ();

my @TUMOR_array_length = ();
my @TUMOR_array = ();

my $k = 0;
my $i = 0;
my $j = 0;

## my @BaseCounts = ();
my @DEPTH = ();

my @FS = ();
my @SOR = ();

###################### setting the variables for FORMAT :
###################### NORMAL and TUMOR samples .

my @NORMAL_GT = ();
my @NORMAL_AD = ();
my @NORMAL_AD_REF = ();
my @NORMAL_AD_ALT = ();
my @NORMAL_AF = ();

my @TUMOR_GT = ();
my @TUMOR_AD = ();
my @TUMOR_AD_REF = ();
my @TUMOR_AD_ALT = ();
my @TUMOR_AF = ();


################################
################################

# Reading the VCF file

while(defined($line = <INPUT>))
{
	chomp $line;
        $lines[$i]=$line;

if ($line =~m/#/g) { print OUTPUT "$lines[$i]\n"; }

if ($line !~m/#/g)
{
	          @t = split(/\t/, $line);

                  $CHROM[$i] = $t[0];
                  $POS[$i] = $t[1];
                  $ID[$i] = $t[2];
                  $REF[$i] = $t[3];
                  $ALT[$i] = $t[4];
                  $QUAL[$i] = $t[5];
                  $FILTER[$i] = $t[6];
                  $INFO[$i] = $t[7];
                  $FORMAT[$i] = $t[8];
                  $NORMAL[$i] = $t[9];
                  $TUMOR[$i] = $t[10];


## we do parse the INFO :

                  @INFO_array = split(/;/, $INFO[$i]);
                  $INFO_array_length[$i] =  @INFO_array;

                  ## print OUTPUT "$lines[$i]\t";

## we do a separate analysis on SNP, INSERTION, DELETION.
## the info on FS and SOR only for SNV.

            if ($INFO[$i] =~m/SNP/g)
            {

                      for ($j=0; $j < $INFO_array_length[$i]; $j++)
                      {
                        if ($INFO_array[$j] =~m/FS=(\d*.\d*)/g)
                         {
                             $FS[$i] = $1 ;
                             ## print OUTPUT "$lines[$i]\t";
                             ## print OUTPUT "$INFO_array[$j] \t $FS[$i] \t";
                         }
                      }

                      for ($j=0; $j < $INFO_array_length[$i]; $j++)
                      {
                        if ($INFO_array[$j] =~m/SOR=(\d*.\d*)/g)
                          {
                            $SOR[$i] = $1 ;
                            ## print OUTPUT "$lines[$i]\t";
                            ## print OUTPUT "$INFO_array[$j] \t $SOR[$i] \t";
                          }
                      }

                     $number_SNP_unfiltered = $number_SNP_unfiltered + 1;

              }

            if ($INFO[$i] =~m/INSERTION/g)
            {
                     $number_INSERTION_unfiltered = $number_INSERTION_unfiltered + 1;
            }


            if ($INFO[$i] =~m/DELETION/g)
            {
                     $number_DELETION_unfiltered = $number_DELETION_unfiltered + 1;
            }

################################################################################
##### we parse the other fields : FORMAT, NORMAL and TUMOR.
################################################################################

                   ## print OUTPUT "$lines[$i]\t";

                   @FORMAT_array = split(/:/, $FORMAT[$i]);
                   $FORMAT_array_length[$i] =  @FORMAT_array;

                   @NORMAL_array = split(/:/, $NORMAL[$i]);
                   $NORMAL_array_length[$i] =  @NORMAL_array;

                   @TUMOR_array = split(/:/, $TUMOR[$i]);
                   $TUMOR_array_length[$i] =  @TUMOR_array;

###################################################################################
###################################################################################

                   # $NORMAL_GT[$i]  = $NORMAL_array[0];
                   # $NORMAL_AD[$i]  = $NORMAL_array[1];
                   # $NORMAL_AF[$i] = $NORMAL_array[2];

                   # if ($NORMAL_AD[$i] =~m/(\d*),(\d*)/g)
                   # {
                   #    $NORMAL_AD_REF[$i] = $1;
                   #    $NORMAL_AD_ALT[$i] = $2;
                   # }

## print OUTPUT "$NORMAL_GT[$i]\t$NORMAL_AD[$i]\t$NORMAL_AF[$i]\t$NORMAL_AD_REF[$i]\t$NORMAL_AD_ALT[$i]\t";

                   $TUMOR_GT[$i] = $TUMOR_array[0];
                   $TUMOR_AD[$i] = $TUMOR_array[1];
                   $TUMOR_AF[$i] = $TUMOR_array[2];

                   if ($TUMOR_AD[$i] =~m/(\d*),(\d*)/g)
                   {
                       $TUMOR_AD_REF[$i]  = $1;
                       $TUMOR_AD_ALT[$i]  = $2;

                       ## here we are computing the depth based on the SUM : AD_REF + AD_ALT
                       $DEPTH[$i] = $TUMOR_AD_REF[$i] + $TUMOR_AD_ALT[$i];
                   }

## print OUTPUT "$TUMOR_GT[$i]\t$TUMOR_AD[$i]\t$TUMOR_AF[$i]\t$TUMOR_AD_REF[$i]\t$TUMOR_AD_ALT[$i]\t";

## we are applying the filtering criteria, separately on SNP and INDEL :
## if there is an INSERTION or DELETION :

if ($INFO[$i] =~ m/VariantType=INSERTION/)
{
             if ((min($TUMOR_AD_ALT[$i],$TUMOR_AD_REF[$i]) >= $MIN_AD) && ($TUMOR_AF[$i] >= $MIN_AF))
             {
                  print OUTPUT "$lines[$i]\n";
                  $number_INSERTION_filtered = $number_INSERTION_filtered + 1 ;
                  ## print "$INFO[$i]\n";
             }
}

if ($INFO[$i] =~ m/VariantType=DELETION/)
{
             if ((min($TUMOR_AD_ALT[$i],$TUMOR_AD_REF[$i]) >= $MIN_AD) && ($TUMOR_AF[$i] >= $MIN_AF))
             {
                  print OUTPUT "$lines[$i]\n";
                  $number_DELETION_filtered = $number_DELETION_filtered + 1 ;
                  ## print "$INFO[$i]\n";
             }
}

if ($INFO[$i] =~ m/VariantType=SNP/)
{
             if ( (min($TUMOR_AD_ALT[$i],$TUMOR_AD_REF[$i]) >= $MIN_AD) && ($TUMOR_AF[$i] >= $MIN_AF) &&
                                                           ($FS[$i] < $MAX_FS) && ($SOR[$i] < $MAX_SOR) )
             {
                  print OUTPUT "$lines[$i]\n";
                  $number_SNP_filtered = $number_SNP_filtered + 1;
                  ## print "$INFO[$i]\n";
                  ## print OUTPUT "$FS[$i] \t";
                  ## print OUTPUT "$SOR[$i] \n";
             }
}


}

}

$number_VARIANTS_unfiltered = $number_SNP_unfiltered + $number_INSERTION_unfiltered + $number_DELETION_unfiltered ;

print OUTPUT2 "\n";
print OUTPUT2 "TOTAL number of un-filtered VARIANTS is :\t $number_VARIANTS_unfiltered \n";
print OUTPUT2 "TOTAL number of un-filtered SNPs is :\t $number_SNP_unfiltered \n";
print OUTPUT2 "TOTAL number of un-filtered INSERTIONS is :\t $number_INSERTION_unfiltered \n";
print OUTPUT2 "TOTAL number of un-filtered DELETIONS is :\t $number_DELETION_unfiltered \n";

$number_VARIANTS_filtered = $number_SNP_filtered + $number_INSERTION_filtered + $number_DELETION_filtered ;

print OUTPUT2 "\n";
print OUTPUT2 "TOTAL number of filtered VARIANTS is :\t $number_VARIANTS_filtered \n";
print OUTPUT2 "TOTAL number of filtered SNPs is :\t $number_SNP_filtered \n";
print OUTPUT2 "TOTAL number of filtered INSERTIONS is :\t $number_INSERTION_filtered \n";
print OUTPUT2 "TOTAL number of filtered DELETIONS is :\t $number_DELETION_filtered \n";

close (INPUT);
close (OUTPUT);
exit;
