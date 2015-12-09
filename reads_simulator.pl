#!/usr/bin/perl

$| = 1;
use FindBin '$Bin';
use Switch;
$path = $Bin;

if(@ARGV < 1) {
    print "\nusage: reads_simulator.pl <num_reads> <name> [options]\n\n";
    print "<name> string is used in the names of the output files.\nUse only alphanumeric, underscores and dashes in this string.\n";

    print "\nThis program outputs a fasta file of reads and a .cig file and .bed file representing\nthe truth of where these reads map.\n";
    print "   - bed file is in one-based coords and contains both endpoints of each span.\n";
    print "   - cig file has a cigar string representation of the mapping coordinates, and a more human\n     readable representation of the coordinates.  All coords are one-based.\n";
    print "Also output are: a file of indels, a file of substitutions, a file of splice forms, and a file of junctions \n";
    print "and gaps crossed by the simulated reads.\n";
    print "\n   options:\n";
    print "      -readlength n   : Set readlength to n>0 bases (default n = 100).\n";
    print "      -fraglength n,k,m : Set minimum fragment length to n bases and the max to m bases and median to k bases.  m > k > n and n must be at least as long as the readlength.\n";
    print "      -numgenes n     : Choose n>=1 genes at random from a master pool of gene models (default n = 30000).\n                        This only works with the four config files that have no stem (see below for more info).\n";
    print "      -strandspecific : Generate strand specific data\n";
    print "      -error x        : Set the error rate that any given base is sequenced wrong to 0<=x<=1 (default x = 0.005).\n";
    print "      -subfreq x      : Set substitution rate to 0<=x<1 (default x = 0.001).\n";
    print "      -indelfreq x    : Set indel rate to 0<=x<1 (default x = 0.0005).\n";
    print "      -intronfreq x   : Set intron rate to 0<=x<1 (default is determined from the\n";
    print "                        featureqantification file).\n";
    print "      -tlen n         : Set the length of the low quality tail to n bases (default n = 10).\n";
    print "      -tpercent x     : Set the percent of tails that are low quality to 0<=x<=1\n                        (default x = 0).\n";
    print "      -tqual x        : Set quality of the low quality tail to 0<=x<=1 (default x = 0.8).\n";
    print "      -nalt n         : Set the number of novel splice forms per gene to n>1 (default n = 2).\n";
    print "      -palt x         : Set the percentage of signal coming from novel splice\n                        forms to 0<=x<=1 (default x = 0.2).\n";
    print "      -sn             : Add sequence number to the first column of the bed file.\n";
    print "      -configstem x   : The stem x will be added to the four default filenames for the four\n                        required files.  Use this if you have made your own config files (see\n                        below about config files and custom config files).\n";
    print "                        E.g. \"simulator_config_geneinfo\" becomes \"simulator_config_geneinfo_x\".\n";
    print "      -cntstart n     : Start the read counter at n (default n = 1).\n";
    print "      -outdir x       : x is a path to a directory to write to.  Default is the current directory.\n";
    print "      -mastercfgdir x : x is a path to a directory where the master config files are.  Default is the\n                        current directory.\n";
    print "      -customcfgdir x : If you are using -configstem option, then x is a path to a directory where the\n";
    print "                        custom config files are.  Default is the directory specified by -outdir, which.\n";
    print "                        itself defaults to the current directory.\n";
    print "      -usesubs x      : x is a file of substitutions in the format output by this program in case you want\n";
    print "                        to resuse them for another run with the same gene models.\n";
    print "      -useindels x    : x is a file of indels in the format output by this program in case you want\n";
    print "                        to resuse them for another run with the same gene models.\n";
    print "      -usealts x      : x is a file of transcripts in the format output by this program in case you want\n";
    print "                        to resuse them for another run with the same gene models.  This must be used in\n                        conjunction with -configstem for the same four files used the first time (those\n                        will be the files ending in 'temp' if you used the master pool the first time).\n";
    print "\n";
    print "This program depends on four (master config) files:\n";
    print "  1) simulator_config_geneinfo\n";
    print "  2) simulator_config_geneseq\n";
    print "  3) simulator_config_intronseq\n";
    print "  4) simulator_config_featurequantifications\n\n";
    print "With these files the program chooses -numgenes of them at random.\nIf you create custom config files with a suffix and use the -configstem mode\nit then uses all genes in the file.\n";
    print "To create custom config files for a subset of genes in the master config, use the script:\n";
    print "   - make_config_files_for_subset_of_gene_ids.pl\n";
    print "Run it with no parameters for the usage\n";
    print "To use the custom config files with this program put them in the same directory as the script, or\n";
    print "in the directory specified by -outdir (not -mastercfgdir which specifies the master config files)\n";
    print "and use the option -configstem\n";
    print "\n";
    exit(0);
}

$name = $ARGV[1];

if($name =~ /^(-|_)/) {
    print STDERR "\nError: invalid name '$name' - name cannot start with a dash or an underscore.\n\n";
    exit(0);
}
if($name =~ /[^a-zA-Z0-9_\-\.]/) {
    print STDERR "\nError: name must be only letters, numbers, dashes, underscores and periods.\n\n";
    exit(0);
}
if(!($name =~ /[a-zA-Z0-9_\-\.]/)) {
    print STDERR "\nERROR: invalid name '$name', cannot be empty and must use only letters,\n       numbers, dashes, underscores and periods.\n\n";
    exit(0);
}

$READLENGTH = 100;
$FRAGLENGTH = "";
$substitutionfrequency = .001;
$indelfrequency = .0005;
$base_error = .005;
$low_qual_tail_length = 10;
$percent_of_tails_that_are_low_qual = 0;
$quality_of_low_qual_tail = .8;
$percent_alt_spliceforms = .2;
$num_alt_splice_forms_per_gene = 2;
$NUMGENES = 30000;
$stem = "";
$cntstart = 1;
$outdir = "./";
$mastercfgdir = "./";
$customcfgdir = "";
$intronfrequency = -1;

BEGIN {
    eval {
        require Math::Random;

        # If you'd ordinarily "use Module::Name qw(foo bar baz);", pass                                          
        # the qw(foo bar baz) to import here.                                                                    

        import Math::Random qw(:all);
    };

    # If the eval failed, we don't have the module                                                               
    if ($@) {
        print STDERR "\nOops, you must first install the Math::Random module for this to work...\n\n";
        exit(0);
    }
}

use Math::Random qw(:all);
$num_reads = $ARGV[0];
if(!($num_reads =~ /^\d+$/) || ($num_reads <= 0)) {
    print STDERR "\nError: number of reads must be an integer, not '$ARGV[0]'.\n\n";
    exit(0);
}
$sum_of_gene_counts = 0;
$sum_of_intron_counts = 0;
$genecnt=0;
$exoncount_total=0;
$introncount_total=0;

$usesubs = "false";   # will become true if there is a custom subsitutions file specified
$useindels = "false"; # will become true if there is a custom indels file specified
$usealts = "false"; # will become true if there is a custom transcripts file specified
$seq_num_in_bedfile = "false";
$strandspecific = "false";
$outputfq = "false";
$varcov = "";
$dampen = 1;
$usevarcov = "false";
$usepcd = "false";
for($i=2; $i<@ARGV; $i++) {
    $option_recognized = 0;
    if($ARGV[$i] eq "-sn") {
	$seq_num_in_bedfile = "true";
	$option_recognized = 1;
    }
    if($ARGV[$i] eq "-outputfq") {
	$outputfq = "true";
	$option_recognized = 1;
    }
    if($ARGV[$i] eq "-strandspecific") {
	$strandspecific = "true";
	$option_recognized = 1;
    }
    if($ARGV[$i] eq "-configstem") {
	$i++;
	$stem = $ARGV[$i];
	$option_recognized = 1;
	if($stem =~ /^(-|_)/) {
	    print STDERR "\nError: -configstem cannot start with a dash or an underscore.\n\n";
	    exit(0);
	}
	if($stem =~ /[^a-zA-Z0-9._-]/) {
	    print STDERR "\nError: -configstem must be only letters, numbers, dashes, periods and underscores.\n\n";
	    exit(0);
	}

	$simulator_config_geneinfo = "simulator_config_geneinfo_$stem";
	$simulator_config_featurequantifications = "simulator_config_featurequantifications_$stem";
	$simulator_config_geneseq = "simulator_config_geneseq_$stem";
	$simulator_config_intronseq = "simulator_config_intronseq_$stem";
    }
    if($ARGV[$i] eq "-tqual") {
	$i++;
	$quality_of_low_qual_tail = $ARGV[$i];
	$option_recognized = 1;
	if(!($quality_of_low_qual_tail =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -tqual has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	$quality_of_low_qual_tail = $quality_of_low_qual_tail + 0;
	if($quality_of_low_qual_tail > 1 || $quality_of_low_qual_tail < 0) {
	    print STDERR "\nError: -tqual has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-tpercent") {
	$i++;
	$percent_of_tails_that_are_low_qual = $ARGV[$i];
	$option_recognized = 1;
	if(!($percent_of_tails_that_are_low_qual =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -tpercent has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	$percent_of_tails_that_are_low_qual = $percent_of_tails_that_are_low_qual + 0;
	if($percent_of_tails_that_are_low_qual > 1 || $percent_of_tails_that_are_low_qual < 0) {
	    print STDERR "\nError: -tpercent has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-tlen") {
	$i++;
	$low_qual_tail_length = $ARGV[$i];
	$option_recognized = 1;
	if(!($low_qual_tail_length =~ /^\d+$/)) {
	    print STDERR "\nError: -teln must be a positive integer.\n\n";
	    exit(0);
	}
	if($low_qual_tail_length < 1) {
	    print STDERR "\nError: -teln must be a positive integer.\nIf you want there to be no tail error don't set -tlen\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-error") {
	$i++;
	$base_error = $ARGV[$i];
	$option_recognized = 1;
	if(!($base_error =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -error has to be non-negative and less than one\n\n";
	    exit(0);
	}
	$base_error = $base_error + 0;
	if($base_error >= 1 || $base_error < 0) {
	    print STDERR "\nError: -error has to be non-negative and less than one\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-nalt") {
	$i++;
	$num_alt_splice_forms_per_gene = $ARGV[$i];
	$option_recognized = 1;
	if(!($num_alt_splice_forms_per_gene =~ /^\d+$/)) {
	    print STDERR "\nError: -nalt has to be a postive integer.\nIf you want there to be no alternate splice forms don't set -nalt, set -palt 0\n\n";
	    exit(0);
	}
	if($num_alt_splice_forms_per_gene < 1) {
	    print STDERR "\nError: -nalt has to be a postive integer.\nIf you want there to be no alternate splice forms don't set -nalt, set -palt 0\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-numgenes") {
	$i++;
	$NUMGENES = $ARGV[$i];
	$option_recognized = 1;
	if(!($NUMGENES =~ /^\d+$/)) {
	    print STDERR "\nError: -numgenes has to be a postive integer.\n\n";
	    exit(0);
	}
	if($NUMGENES < 1) {
	    print STDERR "\nError: -numgenes has to be a postive integer.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-cntstart") {
	$i++;
	$cntstart = $ARGV[$i];
	$option_recognized = 1;
	if(!($cntstart =~ /^\d+$/)) {
	    print STDERR "\nError: -cntstart has to be a postive integer.\n\n";
	    exit(0);
	}
	if($cntstart < 1) {
	    print STDERR "\nError: -cntstart has to be a postive integer.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-outdir") {
	$i++;
	$outdir = $ARGV[$i];
	$outdir =~ s!/$!!;
	$outdir = $outdir . "/";
	$option_recognized = 1;
	if(!(-e $outdir)) {
	    print STDERR "\nError: cannot open the directory '$outdir' specified by the -outdir option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-varcov") {
	$i++;
	$varcov = $ARGV[$i];
	$option_recognized = 1;
	if(!(-e $varcov)) {
	    print STDERR "\nError: cannot open the file '$varcov' specified by the -varcov option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-dampen") {
	$i++;
	$dampen = $ARGV[$i];
	$option_recognized = 1;
	if(!($dampen =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -dampen has to be a number between 0 and 1 (inclusive).\n\n";
	    exit(0);
	}
	$dampen = $dampen + 0;
	if($dampen > 1 || $dampen < 0) {
	    print STDERR "\nError: -dampen has to be strictly between 0 and 1 (inclusive).\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-mastercfgdir") {
	$i++;
	$mastercfgdir = $ARGV[$i];
	$mastercfgdir =~ s!/$!!;
	$mastercfgdir = $mastercfgdir . "/";
	$option_recognized = 1;
	if(!(-e $mastercfgdir)) {
	    print STDERR "\nError: cannot open the directory '$mastercfgdir' specified by the -mastercfgdir option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-usepcd") {
	$option_recognized = 1;
	$usepcd = "true";
    }
    if($ARGV[$i] eq "-usesubs") {
	$i++;
	$subsfile = $ARGV[$i];
	$option_recognized = 1;
	$usesubs = "true";
	if(!(-e $subsfile)) {
	    print STDERR "\nError: cannot open the file '$subsfile' specified by the -usesubs option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-useindels") {
	$i++;
	$indelsfile = $ARGV[$i];
	$useindels = "true";
	$option_recognized = 1;
	if(!(-e $indelsfile)) {
	    print STDERR "\nError: cannot open the file '$indelsfile' specified by the -useindels option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-usealts") {
	$i++;
	$altsfile = $ARGV[$i];
	$usealts = "true";
	$option_recognized = 1;
	if(!(-e $altsfile)) {
	    print STDERR "\nError: cannot open the file '$altsfile' specified by the -usealts option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-usevarcov") {
	$i++;
	$varcovfile = $ARGV[$i];
	$option_recognized = 1;
	$usevarcov = "true";
	if(!(-e $varcovfile)) {
	    print STDERR "\nWarning: Cannot open the file '$varcovfile' specified by the -usevarcov option,\ngoing to assume you meant not to use this option.\n\n";
	    $usevarcov = "false";
	    $usepcd = "false";
	}
    }
    if($ARGV[$i] eq "-customcfgdir") {
	$i++;
	$customcfgdir = $ARGV[$i];
	$customcfgdir =~ s!/$!!;
	$customcfgdir = $customcfgdir . "/";
	$option_recognized = 1;
	if(!(-e $customcfgdir)) {
	    print STDERR "\nError: cannot open the directory '$customcfgdir' specified by the -customcfgdir option.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-palt") {
	$i++;
	$percent_alt_spliceforms = $ARGV[$i];
	$option_recognized = 1;
	if(!($percent_alt_spliceforms =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -palt has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	$percent_alt_spliceforms = $percent_alt_spliceforms + 0;
	if($percent_alt_spliceforms > 1 || $percent_alt_spliceforms < 0) {
	    print STDERR "\nError: -palt has to be non-negative and no more than one.\n\n";
	    exit(0);
	}
	if($percent_alt_spliceforms == 0) {
	    $num_alt_splice_forms_per_gene = 0;
	}
    }
    if($ARGV[$i] eq "-indelfreq") {
	$i++;
	$indelfrequency = $ARGV[$i];
	$option_recognized = 1;
	if(!($indelfrequency =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -indelfreq has to be strictly between 0 and 1.\n\n";
	    exit(0);
	}
	$indelfrequency = $indelfrequency + 0;
	if($indelfrequency >= 1 || $indelfrequency < 0) {
	    print STDERR "\nError: -indelfreq has to be strictly between 0 and 1.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-intronfreq") {
	$i++;
	$intronfrequency = $ARGV[$i];
	$option_recognized = 1;
	if(!($intronfrequency =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -intronfreq has to be between 0 and 1.\n\n";
	    exit(0);
	}
	$intronfrequency = $intronfrequency + 0;
	if($intronfrequency >= 1 || $intronfrequency < 0) {
	    print STDERR "\nError: -intronfreq has to be between 0 and 1.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-subfreq") {
	$i++;
	$substitutionfrequency = $ARGV[$i];
	$option_recognized = 1;
	if(!($substitutionfrequency =~ /^\d*\.?\d*$/)) {
	    print STDERR "\nError: -subfreq has to be less than one and non-negative.\n\n";
	    exit(0);
	}
	$substitutionfrequency = $substitutionfrequency + 0;
	if($substitutionfrequency >= 1 || $substitutionfrequency < 0) {
	    print STDERR "\nError: -subfreq has to be less than one and non-negative.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-readlength") {
	$i++;
	$READLENGTH = $ARGV[$i];
	$option_recognized = 1;
	if(!($READLENGTH =~ /^\d+$/)) {
	    print STDERR "\nError: -readlength has to be a postive integer.\n\n";
	    exit(0);
	}
	$READLENGTH = $READLENGTH + 0;
	if($READLENGTH < 1) {
	    print STDERR "\nError: -readlength has to be a postive integer.\n\n";
	    exit(0);
	}
    }
    if($ARGV[$i] eq "-fraglength") {
	$i++;
	$FRAGLENGTH = $ARGV[$i];
	$option_recognized = 1;
	if(!($FRAGLENGTH =~ /^\d+,\d+,\d+$/)) {
	    print STDERR "\nError: -fraglength has to be a pair of postive integers separated by a comma, with no whitespace.\n\n";
	    exit(0);
	}
	$FRAGLENGTH =~ /^(\d+),(\d+),(\d+)$/;
	$minfraglength = $1;
	$medfraglength = $2;
	$maxfraglength = $3;
    }
    if($option_recognized == 0) {
	print STDERR "\nError: option $ARGV[$i] not recognized.\n\n";
	exit(0);
    }
}

if($usevarcov eq "true" && (!$varcov =~ /\S/)) {
    print STDERR "\nError: If you use the -usevarcov option which specifies the file mapping sim genes to ivt genes,\nyou still must also specify the original varcov file itself using the -varcov option.\n\n";
    exit(0);
}
if($usepcd eq "true" && $usevarcov eq "false") {
    print STDERR "\nError: If you use the -usepcd option you must also be using the -usevarcov option.\n\n";
    exit(0);
}

#if($minfraglength =~ /\S/ && $minfraglength < $READLENGTH) {
#    print STDERR "\nError: minimum fraglength has to be at least as big as the read length.\n\n";
#    exit(0);
#}

if(!($minfraglength =~ /\S/)) {
    $minfraglength = $READLENGTH;
    $medfraglength = $READLENGTH * 2;
    $maxfraglength = $READLENGTH * 5;
}

if($usesubs eq "false" && $useindels eq "true") {
    die "Error: specify either both -usesubs and -useindels or neither.\n\n";
}
if($usesubs eq "true" && $useindels eq "false") {
    die "Error: specify either both -usesubs and -useindels or neither.\n\n";
}

if($customcfgdir =~ /\S/ && !($stem =~ /\S/)) {
    die "\nError: if you specify -customcfgdir you must also specify -stem.\n\n";
}

if($customcfgdir =~ /\S/) {
    $simulator_config_geneinfo = $customcfgdir . $simulator_config_geneinfo;
    $simulator_config_featurequantifications = $customcfgdir . $simulator_config_featurequantifications;
    $simulator_config_geneseq = $customcfgdir . $simulator_config_geneseq;
    $simulator_config_intronseq = $customcfgdir . $simulator_config_intronseq;
} else {
    $simulator_config_geneinfo = $outdir . $simulator_config_geneinfo;
    $simulator_config_featurequantifications = $outdir . $simulator_config_featurequantifications;
    $simulator_config_geneseq = $outdir . $simulator_config_geneseq;
    $simulator_config_intronseq = $outdir . $simulator_config_intronseq;
}

if($READLENGTH <= $low_qual_tail_length && $percent_of_tails_that_are_low_qual > 0) {
    print STDERR "\nERROR: low quality tail length must be less than the readlength.\n\n";
    exit(0);
}
if(!($stem =~ /\S/)) {
    # Here construct random set of $NUMGENES genes and call the 'make_config_files_for_subset_of_gene_ids'
    # script on the master-list files
    $F = $mastercfgdir . 'simulator_config_geneinfo';
    $total_num_genes = `wc -l $F`;
    $total_num_genes =~ /\s*(\d+)/;
    $total_num_genes = $1;
    $genes_temp_filename = $outdir . "genes_temp";
    open(OUTFILE, ">$genes_temp_filename") or die "\nError: cannot open file '$genes_temp_filename' for writing\n\n";
    for($i=0; $i<$NUMGENES; $i++) {
	$R = int(rand($total_num_genes))+1;
	while($Ghash{$R}+0>0) {
	    $R = int(rand($total_num_genes));
	}
	$Ghash{$R}++;
	print OUTFILE "GENE.$R\n";
    }
    close(OUTFILE);
    $x = `perl $path/make_config_files_for_subset_of_gene_ids.pl temp $genes_temp_filename $mastercfgdir $outdir`;

    $simulator_config_geneinfo = $outdir . "simulator_config_geneinfo_temp";
    $simulator_config_featurequantifications = $outdir . "simulator_config_featurequantifications_temp";
    $simulator_config_geneseq = $outdir . "simulator_config_geneseq_temp";
    $simulator_config_intronseq = $outdir . "simulator_config_intronseq_temp";
}

$logfilename = $outdir . "simulated_reads_$name" . ".log";
$bedfilename = $outdir . "simulated_reads_$name" . ".bed";
$fafilename_forward = $outdir . "simulated_reads_$name" . ".forward.fa";
$fafilename_reverse = $outdir . "simulated_reads_$name" . ".reverse.fa";
if($outputfq eq "true") {
    $fafilename_forward =~ s/fa$/fq/;
    $fafilename_reverse =~ s/fa$/fq/;
}
$substitutionsfilename = $outdir . "simulated_reads_substitutions_$name" . ".txt";
$indelsfilename = $outdir . "simulated_reads_indels_$name" . ".txt";
$junctionssfilename = $outdir . "simulated_reads_junctions-crossed_$name" . ".txt";
$reads2genesfilename = $outdir . "simulated_reads2genes_$name" . ".txt";
$individual_forward_and_reverse_alignments_file = $outdir . "simulated_reads_$name" . ".cig";
#$frags_file = $outdir . "simulated_fragments_$name" . ".txt";

open(SIMLOGOUT, ">$logfilename") or die "\nError: cannot open file '$logfilename' for writing\n\n";
open(SIMBEDOUT, ">$bedfilename") or die "\nError: cannot open file '$bedfilename' for writing\n\n";
open(SIMFAOUTF, ">$fafilename_forward") or die "\nError: cannot open file '$fafilename_forward' for writing\n\n";
open(SIMFAOUTR, ">$fafilename_reverse") or die "\nError: cannot open file '$fafilename_reverse' for writing\n\n";
open(SIMSUBSOUT, ">$substitutionsfilename") or die "\nError: cannot open file '$substitutionsfilename' for writing\n\n";
open(SIMINDELSOUT, ">$indelsfilename") or die "\nError: cannot open file '$indelsfilename' for writing\n\n";
open(SIMJUNCTIONSOUT, ">$junctionssfilename") or die "\nError: cannot open file '$substitutionsfilename' for writing\n\n";
open(SIMREAD2GENE, ">$reads2genesfilename") or die "\nError: cannot open file '$reads2genesfilename' for writing\n\n";
open(SIMCIGOUT, ">$individual_forward_and_reverse_alignments_file") or die "\nError: cannot open file '$individual_forward_and_reverse_alignments_file' for writing\n\n";
#open(FRAGSEQOUT, ">$frags_file") or die "\nError: cannot open file '$frags_file' for writing\n\n";

$date = `date`;
chomp($date);
print SIMLOGOUT "Simulator run: '$name'\n";
print SIMLOGOUT "started: $date\n";
$f = format_large_int($num_reads);
print SIMLOGOUT "num reads: $f\n";
print SIMLOGOUT "readlength: $READLENGTH\n";
print SIMLOGOUT "substitution frequency: $substitutionfrequency\n";
print SIMLOGOUT "indel frequency: $indelfrequency\n";
print SIMLOGOUT "base error: $base_error\n";
if($intronfrequency >= 0) {
    print SIMLOGOUT "Intron frequency: $intronfrequency\n";
} else {
    print SIMLOGOUT "Intron frequency: taken from featurequantification config file\n";
}

print SIMLOGOUT "low quality tail length: $low_qual_tail_length\n";
print SIMLOGOUT "percent of tails that are low quality: $percent_of_tails_that_are_low_qual\n";
print SIMLOGOUT "quality of low qulaity tails: $quality_of_low_qual_tail\n";
print SIMLOGOUT "percent of alt splice forms: $percent_alt_spliceforms\n";
print SIMLOGOUT "number of alt splice forms per gene: $num_alt_splice_forms_per_gene\n";
print SIMLOGOUT "fraglength: $minfraglength, $medfraglength, $maxfraglength\n";
if($stem =~ /\S/) {
    print SIMLOGOUT "stem: $stem\n";
} else {
    $f = format_large_int($NUMGENES);
    print SIMLOGOUT "num genes used to simulate with: $f\n";
}
close(SIMLOGOUT);
open(SIMLOGOUT, ">>$logfilename");

print STDERR "reading the gene info file\n";
open(INFILE, $simulator_config_geneinfo) or die "\nError: cannot open file '$simulator_config_geneinfo' for reading\n\n";
$total_transcriptome_length = 0;
while($line = <INFILE>) {
    chomp($line);
    @a = split(/\t/,$line);
    @b = split(/::::/,$a[7]);
    for($i=0; $i<@b; $i++) {
	$b[$i] =~ s/\(.*\)//;
	$a[5] =~ s/,\s*$//;
	$a[6] =~ s/,\s*$//;
	$starts{$b[$i]} = $a[5];
	$ends{$b[$i]} = $a[6];
	@c1 = split(/,/,$a[5]);
	@c2 = split(/,/,$a[6]);
	my $length=0;
	for(my $j=0; $j<@c1; $j++) {
	    $length = $length + $c2[$j]-$c1[$j];
	}
	$gene_length{$b[$i]} = $length;
	$total_transcriptome_length = $total_transcriptome_length + $length;
	$chr{$b[$i]} = $a[0];
	$strand{$b[$i]} = $a[1];
    }
}
close(INFILE);

print STDERR "reading the feature quantification file\n";
open(INFILE, $simulator_config_featurequantifications) or die "\nError: cannot open file '$simulator_config_featurequantifications' for reading\n\n";
while($line = <INFILE>) {
    chomp($line);
    if($line =~ /(exon)/) {
	@a = split(/\t/,$line);
	@a2 = @{$exon2gene{$a[1]}};
	$N = @a2 + 0;
	$exon2gene{$a[1]}[$N] = $geneid;
	$gene2exon{$geneid}[$exoncnt] = $a[1];
	$exon2count{$a[1]} = $a[2];
	$exoncnt++;
    }
    if($line =~ /(intron)/) {
	@a = split(/\t/,$line);
	@a2 = @{$intron2gene{$a[1]}};
	$N = @a2 + 0;
	$intron2gene{$a[1]}[$N] = $geneid;
	$gene2intron{$geneid}[$introncnt] = $a[1];
	$intron2count{$a[1]} = $a[2];
	$introncnt++;
	$gene2introncnt{$geneid}{$a[1]} = $introncnt;
    }
    if($line =~ /--------/) {
	$line = <INFILE>;
	chomp($line);
	$line =~ s/\(.*//;
	$line =~ s/\s*(\+|-)$//;
	$geneid = $line;
	$exoncnt=0;
	$introncnt=0;
	$line = <INFILE>;
	$line = <INFILE>;
	chomp($line);
	@a = split(/\t/,$line);
	$gene_count[$genecnt] = $a[2];
	$gene_count_hash{$geneid} = $a[2];
	$sum_of_gene_counts = $sum_of_gene_counts + $gene_count[$genecnt];
	$genes[$genecnt] = $geneid;
	$genecnt++;
    }
}
close(INFILE);
#$genecnt--;
$numgenes = $genecnt;

if($varcov =~ /\S/) {
    open(IVT, "$varcov");
    $ivt_counter = 0;
    while($line = <IVT>) {
	$line =~ /^>([^ ]*) .*(.)/;
	$IVT_id[$ivt_counter] = $1;
	$IVT_strand[$ivt_counter] = $2;
	$line = <IVT>;
	$IVT_data[$ivt_counter] = $line;
	$IVTid2cnt{$IVT_id[$ivt_counter]} = $ivt_counter;
	$ivt_counter++;
    }
    if($usevarcov eq "true") {
	open(VCF, "$varcovfile");
	$gene_counter = 0;
	while($line = <VCF>) {
	    chomp($line);
	    @a = split(/\t/,$line);
	    $gene2ivt[$gene_counter] = $IVTid2cnt{$a[1]};
	    if($usepcd eq "true") {
		$gene2percentdropoff[$gene_counter] = $a[2];
	    } else {
		$gene2percentdropoff[$gene_counter] = 1;
	    }
	    $gene_counter++;
	}
    } else {
	$ivto = $outdir . "IVT2GENE.txt";
	open(IVTO, ">$ivto");
	for($i=0; $i<$numgenes; $i++) {
	    $rr = int(rand($ivt_counter));
	    my $G = $IVT_id[$rr];
	    $gene2ivt[$i] = $rr;
	    if($dampen < 1) {
		$pcd = rand();
		while($pcd < $dampen) {
		    $pcd = rand();
		}
		$pcd = int($pcd * 1000) / 1000;
	    } else {
		$pcd = 1;
	    }
	    $gene2percentdropoff[$i] = $pcd;
	    $j=$i+1;
	    print IVTO "GENE.$j\t$G\t$pcd\n";
	}
	close(IVTO);
    }
}

$f = format_large_int(int($numgenes));
print "$f genes total\n";
$f = format_large_int(int($sum_of_gene_counts));
print "sum_of_gene_counts = $f\n";
print SIMLOGOUT "sum of gene counts: $f\n";
if($varcov =~ /\S/) {
    print SIMLOGOUT "varcov file used: $varcov\n";
}
if($dampen =~ /\S/) {
    print SIMLOGOUT "dampen: $dampen\n";
}
close(SIMLOGOUT);
open(SIMLOGOUT, ">>$logfilename");

# making alternate (unknown) splice forms
if($percent_alt_spliceforms > 0) {
    if($usealts eq 'false') {
	print STDERR "creating the alternate (unknown) splice forms\n";
    }
}
for($i=0; $i<$numgenes; $i++) {
    # $i is the original gene id
    # gene $genecnt is the same gene as $i but will be modified to be an alternate splice form
    $geneid = $genes[$i];
    @a = @{$gene2exon{$geneid}};
    $exoncnt = @a;
    $original_gene_count = $gene_count[$i];
    # the following adjusts the original count because now it will be split between
    # the normal and alt spice form
    $gene_count[$i] = $original_gene_count * (1-$percent_alt_spliceforms);
    $genes2[$i] = $geneid;
    for($j=0; $j<$num_alt_splice_forms_per_gene; $j++) {
	$geneid_x = $geneid . "_$j";
	$genes2[$genecnt+$j*$numgenes] = $geneid_x;
	$gene_count[$genecnt+$j*$numgenes] = $original_gene_count * $percent_alt_spliceforms/$num_alt_splice_forms_per_gene;
	$exoncnt_x = 0;
	$str = $starts{$geneid};
	$str =~ s/^\s*,\s*//;
	$str =~ s/\s*,\s*$//;
	@S = split(/,/,$str);
	$str = $ends{$geneid};
	$str =~ s/^\s*,\s*//;
	$str =~ s/\s*,\s*$//;
	@E = split(/,/,$str);
	$FLAG1 = 0;
	$FLAG2 = 0;
	while($FLAG1 == 0 || $FLAG2 == 0) {
	    $FLAG1 = 0;
	    $FLAG2 = 0;
	    $starts_new = "";
	    $ends_new = "";
	    $exoncnt_x = 0;
	    delete($gene2exon{$geneid_x});

	    for($k=0; $k<$exoncnt; $k++) {
		$flip = int(rand(2));
		if($flip == 1) {
		    if($strand{$geneid} eq "+") {
			$starts_new = $starts_new . "$S[$k],";
			$ends_new = $ends_new . "$E[$k],";
		    }
		    if($strand{$geneid} eq "-") {
			$starts_new = "$S[$numexons-$k-1]," . $starts_new;
			$ends_new = "$E[$numexons-$k-1]," . $ends_new;
		    }
		    $gene2exon{$geneid_x}[$exoncnt_x] = $gene2exon{$geneid}[$k];
		    $exoncnt_x++;
		    $FLAG1 = 1; # This makes sure we kept at least one exon
		}
		else {
		    $FLAG2 = 1; # This makes sure we skipped at least one exon
		}
	    }
	    if($exoncnt == 1) {
		$FLAG2 = 1;
	    }
	}
	$starts{$geneid_x} = $starts_new;
	$ends{$geneid_x} = $ends_new;
	$strand{$geneid_x} = $strand{$geneid};
	$chr{$geneid_x} = $chr{$geneid};
    }
    $genecnt++;
}
@genes = @genes2;
$genecnt = @genes;

if($usealts eq "true") {
    print STDERR "reading the user specified alternate (unknown) splice forms\n";
    open(INFILE, $altsfile) or die "\nError: cannot open the file '$altsfile' for reading.\n\n";
    $flag = 0;
    $line = <INFILE>;
    chomp($line);
    while($flag == 0) {
	if($line =~ /---------/) {
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ /genes\[(\d+)\] = (.*)/;
	    $i = $1;
	    $geneid = $2;
	    $genes[$i] = $geneid;
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/.*= //;
	    $starts{$geneid} = $line;
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/.*= //;
	    $ends{$geneid} = $line;
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/.*= //;
	    $strand{$geneid} = $line;
	    $line = <INFILE>;
	    chomp($line);
	    $line =~ s/.*= //;
	    $chr{$geneid} = $line;
	    $line = <INFILE>;
	    chomp($line);
	    $j=0;
	    delete($gene2exon{$geneid});
	    until($line =~ /---------/ || $line eq '') {
		$gene2exon{$geneid}[$j] = $line;
		$j++;
		$line = <INFILE>;
		chomp($line);
	    }
	}
	if($line eq '') {
	    $flag = 1;
	}
    }
}

$alttranscriptsfilename = $outdir . "simulated_reads_transcripts_$name" . ".txt";
open(ALTTRANSCRIPTS, ">$alttranscriptsfilename") or die "\nError: cannot open file '$alttranscriptsfilename' for writing.\n\n";
for($i=0; $i<@genes; $i++) {
    $geneid = $genes[$i];
    @a = @{$gene2exon{$geneid}};
    $exoncnt = @a;
    $str1 = $starts{$geneid};
    $str2 = $ends{$geneid};
    $str3 = $strand{$geneid};
    $str4 = $chr{$geneid};
    print ALTTRANSCRIPTS "---------\ngenes[$i] = $geneid\n";
    print ALTTRANSCRIPTS "starts = $str1\n";
    print ALTTRANSCRIPTS "ends = $str2\n";
    print ALTTRANSCRIPTS "strand = $str3\n";
    print ALTTRANSCRIPTS "chr = $str4\n";
    for($j=0; $j<$exoncnt; $j++) {
	print ALTTRANSCRIPTS "$a[$j]\n";
    }
}
close(ALTTRANSCRIPTS);

if($sum_of_gene_counts <= 0) {
    print STDERR "\nERROR: None of the genes have expression level above zero,\ncheck your feature quantification config file.\n\n";
    print SIMLOGOUT "\nERROR: None of the genes have expression level above zero,\ncheck your feature quantification config file.\n\n";
    exit(0);
}

for($i=0; $i<$genecnt; $i++) {
    $gene_density[$i] = $gene_count[$i] / $sum_of_gene_counts;
    $gene_distribution[$i] = $gene_density[$i];
}

$dist_file =  $outdir . "simulated_reads_gene_distribution_" . $name . ".txt";
open(DIST, ">$dist_file");
print DIST "gene_distribution[0] = $gene_distribution[0]\n";
for($i=1; $i<$genecnt; $i++) {
    $gene_distribution[$i] = $gene_distribution[$i] + $gene_distribution[$i-1];
    print DIST "$genes[$i] = $gene_distribution[$i]\n";
}
close(DIST);

$introncount_total=0;
$sum_of_intron_counts = 0;
foreach $intron (keys %intron2gene) {
    $introns[$introncount_total] = $intron;
    $intron_count[$introncount_total] = $intron2count{$intron};
    $sum_of_intron_counts = $sum_of_intron_counts + $intron_count[$introncount_total];
    $introncount_total++;
}

$exoncount_total=0;
$sum_of_exon_counts = 0;
foreach $exon (keys %exon2gene) {
    $exons[$exoncount_total] = $exon;
    $exon_count[$exoncount_total] = $exon2count{$exon};
    $sum_of_exon_counts = $sum_of_exon_counts + $exon_count[$exoncount_total];
    $exoncount_total++;
}

$f = format_large_int(int($sum_of_intron_counts));
print "sum of intron counts = $f\n";
print SIMLOGOUT "sum of intron counts = $f\n";
$f = format_large_int(int($sum_of_exon_counts));
print "sum of exon counts = $f\n";
print SIMLOGOUT "sum of intron counts = $f\n";
# The following formula is based on the fact that introns are not treated as isolated but are padded with
# the rest of the exons.
print "sum_of_exon_counts = $sum_of_exon_counts\n";
$intron_freq = (2 * $sum_of_intron_counts) / ($sum_of_exon_counts + 2 * $sum_of_intron_counts);
$iftemp = $sum_of_intron_counts / ($sum_of_exon_counts + $sum_of_intron_counts);

if($intronfrequency > -1) {
    $intron_freq = $intronfrequency;
} else {
    print "intron frequency: $iftemp\n";
}
print "padded intron freq = $intron_freq\n";
print SIMLOGOUT "intron frequency: $iftemp\n";
print SIMLOGOUT "padded intron frequency: $intron_freq\n";
close(SIMLOGOUT);
open(SIMLOGOUT, ">>$logfilename");

for($i=0; $i<$introncount_total; $i++) {
    if($sum_of_intron_counts>0) {
	$intron_density[$i] = $intron_count[$i] / $sum_of_intron_counts;
    } else {
	$intron_density[$i] = 0;
    }
    $intron_distribution[$i] = $intron_density[$i];
}
for($i=1; $i<$introncount_total; $i++) {
    $intron_distribution[$i] = $intron_distribution[$i] + $intron_distribution[$i-1];
}

#for($i=0; $i<$introncount_total; $i++) {
#    print "intron_distribution[$i] = $intron_distribution[$i]\t$introns[$i]\t$intron_count[$i]\n";
#}

open(INFILE, $simulator_config_geneseq) or die "\nError: cannot open file '$simulator_config_geneseq' for reading.\n\n";

 # this file has seqs on one line and minus strand seqs revserse complemented
 # header line looks like this: >NM_175370:chr1:58714964-58752833_-
$CNT=0;
print "\nReading gene sequences\n";
while($line = <INFILE>) {
    $CNT++;
    if($CNT % 10000 == 0) {
	print "$CNT\n";
    }
    chomp($line);
    $line =~ /:(.*):\d+-\d+/;
    $chr = $1;
    $line =~ s/:.*//;
    $line =~ s/>//;
    $gene_id = $line;
    $SEQ = <INFILE>;
    chomp($SEQ);
    $seq{$gene_id} = $SEQ;
    @S = split(/,/, $starts{$gene_id});
    @E = split(/,/, $ends{$gene_id});
    $offset = 0;
    for($exon=0; $exon<@S; $exon++) {
	$S[$exon]++;
	$len = $E[$exon] - $S[$exon] + 1;
	$EXONSEQ = substr($SEQ,$offset,$len);
	$exonname = $chr . ":" . $S[$exon] . "-" . $E[$exon];
	$exonseq{$exonname} = $EXONSEQ;
	$offset = $offset + $len;
    }
}
close(INFILE);

open(INFILE, $simulator_config_intronseq) or die "\nError: cannot open file '$simulator_config_intronseq' for reading.\n\n";

print "\nReading intron sequences\n";
$flag = 0;
$line = <INFILE>;
$CNT=0;
while($flag == 0) {
    $CNT++;
    if($CNT % 10000 == 0) {
	print "$CNT\n";
    }
    $line =~ /^>(.*)/;
    $intron = $1;
    $seq = "";
    $line = <INFILE>;
    until($line =~ /^>/ || $flag == 1) {
	chomp($line);
	if($line !~ /\S/) {
	    $flag = 1;
	}
	else {
	    if($line !~ /^>/) {
		$seq = $seq . $line;
	    }
	}
	$line = <INFILE>;
    }
    $intronseq{$intron} = $seq;
}
close(INFILE);

# The following puts substitutions and indels into each exon

if($usesubs eq "true") {  # user want to read substitutions from a file
    print STDERR "adding substitutions from the file '$subsfile'\n";
    open(SUBSINFILE, $subsfile);
    while($sub = <SUBSINFILE>) {
	chomp($sub);
	@ssub = split(/\t/,$sub);
	$exon = $ssub[0];
	$exon =~ /^(.*):(\d+)-(\d+)/;
	$start = $2;
	$SEQ = $exonseq{$exon};
	$LOC = $ssub[1] - $start + 1;
	$substitutions_locs{$LOC}++;
	$i = @{$substitutions{$exon}} + 0;
	$substitutions{$exon}[$i] = $LOC;
	$orig = substr($SEQ,$LOC-1,1);
	$B = $ssub[2];
	$B =~ s/.*>//;
	$Z = substr($SEQ,$LOC-1,1,$B);
	$exonseq{$exon} = $SEQ;
	$C = $start + $LOC - 1;
	print SIMSUBSOUT "$exon\t$C\t$Z->$B\n";
    }
    close(SUBSINFILE);
}

if($useindels eq "true") {  # user wants to read indels from a file
    print STDERR "adding indels from the file '$indelsfile'\n";
    open(INDELSINFILE, $indelsfile);
    while($ind = <INDELSINFILE>) {
	chomp($ind);
	@iind = split(/\t/,$ind);
	$exon = $iind[0];
	$SEQ = $exonseq{$exon};
	$LOC = $iind[1];
	$indellength = $iind[2];
	$insert = $iind[3];
	$j = @{$indels{$exon}} + 0;
	$indels{$exon}[$j][0] = $LOC;
	$indels{$exon}[$j][1] = $indellength;
	$indels{$exon}[$j][2] = $insert;
	if($indellength > 0) {  # it's an insertion
	    print SIMINDELSOUT "$exon\t$LOC\t$indellength\t$insert\n";
	    $Z = substr($SEQ,$LOC,0,$insert);
	    $exonseq{$exon} = $SEQ;
	} else {                # it's an deletion
	    $l = -1 * $indellength;
	    $Z = substr($SEQ,$LOC,$l,"");
	    print SIMINDELSOUT "$exon\t$LOC\t$indellength\t$Z\n";
	    $exonseq{$exon} = $SEQ;		    
	}

    }
    close(INDELSINFILE);
}

if($usesubs eq "false" && $useindels eq "false") {  # user wants de novo substitutions and indels, 
                                                      # as opposed to being read from files
    foreach $exon (keys %exon2gene) {
	$try_cnt = 0;
	$exon =~ /^(.*):(\d+)-(\d+)/;
	$chr = $1;
	$start = $2;
	$end = $3;
	$length = $end - $start + 1;
	$num_substitutions = random_binomial(1, $length, $substitutionfrequency);
	undef %substitutions_locs;
	undef @indels_temp;
	$SEQ = $exonseq{$exon};
	for($i=0; $i<$num_substitutions; $i++) {
	    $LOC = int(rand($length)) + 1;
	    while($substitutions_locs{$LOC}+0>0) {
		$LOC = int(rand($length)) + 1;
	    }
	    $substitutions_locs{$LOC}++;
	    $substitutions{$exon}[$i] = $LOC;
	    $orig = substr($SEQ,$LOC-1,1);
	    $B = getrandombase();
	    while($B eq $orig) {
		$B = getrandombase();
	    }
	    $Z = substr($SEQ,$LOC-1,1,$B);
	    $exonseq{$exon} = $SEQ;
	    $C = $start + $LOC - 1;
	    print SIMSUBSOUT "$exon\t$C\t$Z->$B\n";
	}
	
	$num_indels = random_binomial(1, $length, $indelfrequency);
	$flag = 0;
	
	while($flag == 0) {  # the following gets the indel locations and makes sure
	    # indels are at least two bases apart and are in different
	    # places from the substitutions
	    $flag = 1;
	    undef %indels_locs_temp;
	    undef %indels_locs;
	    
	    
	    # first have to choose the indel locations:
	    for($i=0; $i<$num_indels; $i++) {
		$LOC = int(rand($length)+1);
		while(($substitutions_locs{$LOC}+0>0) && ($indels_locs_temp{$LOC}+0>0)) {
		    $LOC = int(rand($length)+1);
		}
		$indels_locs_temp{$LOC}++;
	    }
	    foreach $LOC (keys %indels_locs_temp) {
		if($indels_locs_temp{$LOC} > 0) {
		    $indels_locs{$LOC}++;
		}
	    }
	    foreach $LOC1 (keys %indels_locs) {
		foreach $LOC2 (keys %indels_locs) {
		    $X = $LOC1 - $LOC2;
		    if($X < 2 && $X > -2 && $X != 0) {
			$flag = 0;
		    }
		}	    
	    }	
	    $indelcounter=0;
	    # second have to choose the indel lengths:
	    foreach $LOC (sort {$b<=>$a} keys %indels_locs) {
		$indellength = int(random_exponential(1, 1));
		while($indellength < 1) {
		    $indellength = int(random_exponential(1, 1));
		}
		$flip = int(rand(2));
		if($flip == 1) {                                       # INSERTION
		    $insert = "";
		    for($i=0; $i<$indellength; $i++) {
			$insert = $insert . getrandombase();
		    }
		    $indels_temp[$indelcounter][0] = $LOC;
		    $indels_temp[$indelcounter][1] = $indellength;
		    $indels_temp[$indelcounter][2] = $insert;
		    $indelcounter++;
		}
		else {                                                 # DELETION
		    # have to make sure we don't delete a substitution or part of an insertion or overlap
		    # two deletions, and make sure it doesn't overrun the end of the sequence
		    foreach $LOC2 (keys %substitutions_locs) {
			if($LOC2 <= $LOC + $indellength && $LOC2 > $LOC) {
			    $flag = 0;
			}
		    }
		    foreach $LOC2 (keys %indels_locs) {
			if($LOC2 <= $LOC + $indellength && $LOC2 > $LOC) {
			    $flag = 0;
			}
			if($LOC + $indellength >= $length) {
			    $flag = 0;
			}		    
		    }
		    $indels_temp[$indelcounter][0] = $LOC;
		    $indels_temp[$indelcounter][1] = -1 * $indellength;
		    $indelcounter++;
		}
	    }
	    if($flag == 1) {  # set of indels are in kosher posiitions
		for($j=0;$j<$indelcounter;$j++) {
		    $LOC = $indels_temp[$j][0];
		    $indellength = $indels_temp[$j][1];
		    $insert = $indels_temp[$j][2];
		    $indels{$exon}[$j][0] = $LOC;  # This is the location within the exon, where the
                                                   # first base in the exon is location 1.
		    $indels{$exon}[$j][1] = $indellength;
		    $indels{$exon}[$j][2] = $insert;
		    if($indellength > 0) {
			$Z = substr($SEQ,$LOC,0,$insert);
			$exonseq{$exon} = $SEQ;
			print SIMINDELSOUT "$exon\t$LOC\t$indellength\t$insert\n";
		    }
		    if($indellength < 0) {
			$l = -1 * $indellength;
			$Z = substr($SEQ,$LOC,$l,"");
			$exonseq{$exon} = $SEQ;
			print SIMINDELSOUT "$exon\t$LOC\t$indellength\t$Z\n";
		    }
		}
	    } else {
		$try_cnt++;
		if($try_cnt > 500) { # can't seem to find that many kosherly placed indels
		    # in this exon, going for fewer
		    $num_indels = int($num_indels/2);
		    $try_cnt=0;
		}
	    }	
	}
    }
}
close(SIMINDELSOUT);
close(SIMSUBSOUT);

# Now need to assemble the exons back into genes
print "\nAssembling exons\n";
foreach $geneid (keys %gene2exon) {
    @a = @{$gene2exon{$geneid}};
    $numexons = @a;
    $SEQ = "";
    $indelcounter=0; # counts the number of indels per *gene*
    $offset = 0;
    for($j=0; $j<$numexons; $j++) {
	if($strand{$geneid} eq "+") {
	    $i = $j;
	} else {
	    $i = $numexons - $j - 1;
	}
	$exon = $a[$i];
	$exon =~ /^[^:]+:(\d+)-(\d+)/;
	$exonstart = $1;
	$exonend = $2;
	@a2 = @{$indels{$exon}};
        # need to make hash that associates *genes* to indels in gene specific coords (start of gene = 1)
	for($k=0;$k<@a2;$k++) {
	    $location_in_gene = $offset + $a2[$k][0];
	    $gene2indel_temp{$geneid}[$indelcounter][0] = $location_in_gene;
	    $gene2indel_temp{$geneid}[$indelcounter][1] = $a2[$k][1];  # the length of the indel
	    $indelcounter++;
	}
	$offset = $offset + $exonend - $exonstart + 1;
	$SEQ = $SEQ . $exonseq{$a[$i]};
    }
    $seq{$geneid} = $SEQ;
    $FLAG = 0;
    @a = @{$gene2indel_temp{$geneid}};
    while($FLAG == 0) {
	$FLAG = 1;
	for($i=0; $i<@a-1; $i++) {
	    if($a[$i][0] > $a[$i+1][0]) {
		$temp = $a[$i][0];
		$a[$i][0] = $a[$i+1][0];
		$a[$i+1][0] = $temp;
		$temp = $a[$i][1];
		$a[$i][1] = $a[$i+1][1];
		$a[$i+1][1] = $temp;
		$FLAG = 0;
	    }
	}
    }
    for($i=0; $i<@a; $i++) {
	$gene2indel{$geneid}[$i][0] = $a[$i][0];
	$gene2indel{$geneid}[$i][1] = $a[$i][1];
    }
}

print "\nGenerating reads\n";
$CNT=$cntstart;
$CNT2 = 1;
while( 1 == 1) {
    $R = rand(1);
    if($R < $intron_freq) {
	# pick an intron at random
	$R = rand(1);
	$c = 0;
	while($intron_distribution[$c] < $R && $c<(@intron_distribution-1)) {
	    $c++;
	}
	$INTRON2 = $introns[$c];
	# get a gene at random from the genes containing this intron
	@i2g = @{$intron2gene{$INTRON2}};
	$n = @i2g;
	$totalcount = 0;
	for($i=0; $i<@i2g; $i++) {
	    $totalcount = $totalcount + $gene_count_hash{$i2g[$i]};
	}
	$R2 = int(rand($totalcount))+1;
	$tempcount1 = 0;
	$tempcount2 = 0;
	for($i=0; $i<@i2g; $i++) {
	    $tempcount1 = $tempcount2 + 1;
	    $tempcount2 = $tempcount2 + $gene_count_hash{$i2g[$i]};
	    if($R2 >= $tempcount1 && $R2 <= $tempcount2) {
		$GENE = $i2g[$i];
		$i = @i2g;
	    }
	}
	@g2i = @{$gene2intron{$GENE}};
	$num_introns = @g2i;
	# get the intron number for this intron in this gene
	$intron_num = $gene2introncnt{$GENE}{$INTRON2};

	$SEQ = $seq{$GENE};
	$seqlength = length($SEQ);
	undef %geneWithIntron2indel;
	$SEQ = getpaddedintron($GENE, $intron_num);
	$seqlength = length($SEQ);
	@INDELS = @{$geneWithIntron2indel{$GENE}{$INTRON2}};
# the following fixes the starts/ends so they reflect the retained intron:
	$STARTS2 = $starts{$GENE};
	$STARTS2 =~ s/,\s*$//;
	$STARTS2 =~ s/^\s*,//;
	@S1 = split(/,/, $STARTS2);
	for($ii=0; $ii<@S1; $ii++) {
	    $S1[$ii]++;
	}
	$ENDS2 = $ends{$GENE};
	$ENDS2 =~ s/,\s*$//;
	$ENDS2 =~ s/^\s*,//;
	@E1 = split(/,/, $ENDS2);
	$STARTS = $S1[0];
	$ENDS = "";
	for($p=1; $p<@S1; $p++) {
	    if($strand{$GENE} eq "+") {
		if(!($intron_num == $p)) {
		    $STARTS = $STARTS . ",$S1[$p]";
		    $ENDS = $ENDS . ",$E1[$p-1]";
		}
	    }
	    if($strand{$GENE} eq "-") {
		if(!($intron_num == @S1 - $p)) {
		    $STARTS = $STARTS . ",$S1[$p]";
		    $ENDS = $ENDS . ",$E1[$p-1]";
		}
	    }
	}
	$ENDS = $ENDS . ",$E1[@S1-1]";
	$STARTS =~ s/,\s*$//;
	$STARTS =~ s/^\s*,//;
	$ENDS =~ s/,\s*$//;
	$ENDS =~ s/^\s*,//;
	@S1 = split(/,/, $STARTS);
	$STARTS2 = "";
	for($ii=0; $ii<@S1; $ii++) {
	    $S1[$ii]--;
	    $STARTS2 = $STARTS2 . "$S1[$ii],";
	}
	$STARTS2 =~ s/,\s*$//;
	$STARTS2 =~ s/^\s*,//;
	$STARTS = $STARTS2;

	@S1 = split(/,/, $starts{$GENE});
	@E1 = split(/,/, $ends{$GENE});

	@S1 = split(/,/, $STARTS);
	@E1 = split(/,/, $ENDS);

    }
    else {
	$R = rand(1);
	$c = 0;
	while($gene_distribution[$c] < $R && $c<(@gene_distribution-1)) {
	    $c++;
	}
	$GENE = $genes[$c];
	@INDELS = @{$gene2indel{$GENE}};
	$SEQ = $seq{$GENE};
	$STARTS = $starts{$GENE};
	$ENDS = $ends{$GENE};
    }
    if($strandspecific eq "false") {
	$rev_flip = int(rand(2));
    } else {
	if($strand{$GENE} eq "-") {
	    $rev_flip = 1;
	} else {
	    $rev_flip = 0;
	}
     }

    $return_vector_ref = getreads($SEQ, \@INDELS, $STARTS, $ENDS, $CNT, $chr{$GENE}, $c);
    @return_vector = @{$return_vector_ref};
    $fa = $return_vector[0];
    $bed = $return_vector[1];
    $mergedcoords = $return_vector[2];
    $cigstr1 = $return_vector[3];
    $cigstr2 = $return_vector[4];

    if($cigstr1 =~ /\S/) {
	@C = split(/\t/,$cigstr1);
	$cigstr1 = $C[0];
	$cigstr1 = $cigstr1 . "\t$GENE";
	for($i=1; $i<@C; $i++) {
	    $cigstr1 = $cigstr1 . "\t$C[$i]";
	}
	$cigstr1 = $cigstr1 . "\t";
    }
    if($cigstr2 =~ /\S/) {
	@C = split(/\t/,$cigstr2);
	$cigstr2 = $C[0];
	$cigstr2 = $cigstr2 . "\t$GENE";
	for($i=1; $i<@C; $i++) {
	    $cigstr2 = $cigstr2 . "\t$C[$i]";
	}
	$cigstr2 = $cigstr2 . "\t";
    }

    @falines = split(/\n/,$fa);
    $fatemp = $falines[1];
    $fatemp =~ s/[^N]//g;
    if(length($fatemp) > $readlength / 4) {
	$fa = "none";
    }
    $fatemp = $falines[3];
    $fatemp =~ s/[^N]//g;
    if(length($fatemp) > $readlength / 4) {
	$fa = "none";
    }
    if($fa ne "none") {
	print SIMREAD2GENE "seq.$CNT\t$GENE\n";

	@F = split(/\n/,$fa);
	$return_vector_ref = &add_sequencing_error($F[1], $base_error);
	@return_vector = @{$return_vector_ref};
	$stemp = $return_vector[0];
	$qstring1 = $return_vector[1];

	$return_vector_ref = &add_error_to_tail($stemp, $low_qual_tail_length, $percent_of_tails_that_are_low_qual, $quality_of_low_qual_tail, $qstring1);
	@return_vector = @{$return_vector_ref};
        $stemp2 = $return_vector[0];
	$qstring1 = $return_vector[1];
	$ONE = "$F[0]\n$stemp2\n";
	$S1_temp = $stemp2;

	$return_vector_ref = &add_sequencing_error($F[3], $base_error);
	@return_vector = @{$return_vector_ref};
	$stemp = $return_vector[0];
	$qstring2 = $return_vector[1];

	$return_vector_ref = &add_error_to_tail($stemp, $low_qual_tail_length, $percent_of_tails_that_are_low_qual, $quality_of_low_qual_tail, $qstring2);
	@return_vector = @{$return_vector_ref};
        $stemp2 = $return_vector[0];
	$qstring2 = $return_vector[1];
	$TWO = "$F[2]\n$stemp2\n";
	$S2_temp = $stemp2;

	if($rev_flip == 0) {
	    if($outputfq eq "false") {
		print SIMFAOUTF $ONE;
		print SIMFAOUTR $TWO;
	    } else {
		$ONE_temp = $ONE;
		chomp($ONE_temp);
		$ONE_temp =~ s/\n.*//s;
                $ONE_temp =~ s/>//;
		print SIMFAOUTF "\@ $name $ONE_temp 1\n";
		$ONE_temp = $ONE;
		chomp($ONE_temp);
		$ONE_temp =~ s/.*\n//s;
		print SIMFAOUTF "$ONE_temp\n";
		print SIMFAOUTF "+\n";
		print SIMFAOUTF "$qstring1\n";

		$TWO_temp = $TWO;
                $TWO_temp =~ s/>//;
		chomp($TWO_temp);
		$TWO_temp =~ s/\n.*//s;
		print SIMFAOUTR "\@ $name $TWO_temp 2\n";
		$TWO_temp = $TWO;
		chomp($TWO_temp);
		$TWO_temp =~ s/.*\n//s;
		print SIMFAOUTR "$TWO_temp\n";
		print SIMFAOUTR "+\n";
		print SIMFAOUTR "$qstring2\n";
	    }
	    print SIMCIGOUT "$cigstr1$S1_temp\n";
	    print SIMCIGOUT "$cigstr2$S2_temp\n";
	} else {
	    $ONE =~ s/a/b/;
	    $TWO =~ s/b/a/;
	    if($outputfq eq "false") {
		print SIMFAOUTF $TWO;
		print SIMFAOUTR $ONE;
	    } else {
		$ONE_temp = $ONE;
		chomp($ONE_temp);
		$ONE_temp =~ s/\n.*//s;
                $ONE_temp =~ s/>//;
		print SIMFAOUTR "@ $name $ONE_temp 2\n";
		$ONE_temp = $ONE;
		chomp($ONE_temp);
		$ONE_temp =~ s/.*\n//s;
		print SIMFAOUTR "$ONE_temp\n";
		print SIMFAOUTR "+\n";
		print SIMFAOUTR "$qstring1\n";

		$TWO_temp = $TWO;
		chomp($TWO_temp);
		$TWO_temp =~ s/\n.*//s;
                $TWO_temp =~ s/>//;
		print SIMFAOUTF "@ $name $TWO_temp 1\n";
		$TWO_temp = $TWO;
		chomp($TWO_temp);
		$TWO_temp =~ s/.*\n//s;
		print SIMFAOUTF "$TWO_temp\n";
		print SIMFAOUTF "+\n";
		print SIMFAOUTF "$qstring2\n";
	    }
	    print SIMCIGOUT "$cigstr1$S2_temp\n";
	    print SIMCIGOUT "$cigstr2$S1_temp\n";
	}
	print SIMBEDOUT $bed;
	
	$CNT++;
	$CNT2++;
	if($CNT2 % 50000 == 0) {
	    print "$CNT2 reads done\n";
	}
	if($CNT2 > $num_reads) {
	    close(SIMBEDOUT);
	    close(SIMFAOUTF);
	    close(SIMFAOUTR);
	    close(SIMREAD2GENE);
	    $date = `date`;
	    print SIMLOGOUT "finished at $date\n";
	    close(SIMLOGOUT);

	    $fgl_outfile = $outdir . "fraglenhisto_" . $name . ".txt";
	    open(FGLOUT, ">$fgl_outfile");
	    foreach $key (sort {$a<=>$b} keys %FLhash) {
		print FGLOUT "$key\t$FLhash{$key}\n";
	    }
	    close(FGLOUT);
	    
	    exit(0);
	}
    }
}

close(SIMJUNCTIONSOUT);
#close(FRAGSEQOUT);


sub getreads () {
    ($SEQ, $INDELS_ref, $STARTS, $ENDS, $CNT, $CHR, $genenumber) = @_;

#    The following depends on:
#    1) %seq which maps gene ids to gene sequence
#    2) @INDELS a 2-d array, indel loc is $INDELS[$LOC][0], length is $INDELS[$LOC][1], sorted by $LOC
#    3) $starts and $ends, for the exon start/end positions, zero-based half open
#    4) The current read counter $CNT


    @INDELS = @{$INDELS_ref};

    $seqlength = length($SEQ);
#    if($seqlength < $READLENGTH) {
#	$return_vector[0] = "none";
#	$return_vector[1] = "none";
#	$return_vector[2] = "none";
#	return \@return_vector;
#    }

    $AFlag = 0;
    $ACnt = 0;
    while($AFlag == 0) {
	$ACnt++;
	if($ACnt == 500) {
	    $return_vector[0] = "none";
	    $return_vector[1] = "none";
	    $return_vector[2] = "none";
	    return \@return_vector;
	}
	$AFlag = 1;

	undef @FRAGvect;
    
	for(my $i=1; $i<=$seqlength; $i++) {
	    $FRAGvect[$i] = 0;
	}
	$FRAGvect[0] = 1;
	$FRAGvect[$seqlength+1] = 1;
	$flag = 0;
	
	$numfrags = int($seqlength / $medfraglength);
	if($numfrags <= 2) {
	    $numfrags = 3;
	}
	for($j=0; $j<$numfrags-1; $j++) {
	    $r = int(rand($seqlength));
	    while($FRAGvect[$r] == 1) {
		$r = int(rand($seqlength));
	    }
	    $FRAGvect[$r] = 1;
	}
	$sum = 0;
	$cnt=0;
	$cnt2=0;
	$skipfirstflag = 0;
	
	undef @fraglen;
	undef @fragend;
	for(my $i=1; $i<=$seqlength; $i++) {  # this should skip first and last frags
	    if($FRAGvect[$i] == 0) {
		$cnt++;
	    } else {
		if($cnt >= ($minfraglength - 2) && $cnt <= ($maxfraglength - 3) && $skipfirstflag == 1) {
		    $fraglen[$cnt2] = $cnt;
		    $fragend[$cnt2] = $i;
		    $cnt2++;
		}
		$skipfirstflag = 1;
		$cnt = 0;
	    }
	}
	
	$numfrags = @fraglen;
	if($numfrags == 0) {
	    $AFlag = 0;
	}
	$Flag = 0;
	$RFlag = 0;
	$escape_counter = 0;
	$n = int(random_normal() * 1000 / 3) / 1000;
	if($n < 0) {
	    $n = $n * -1;
	}
	undef %fash;
	while($Flag == 0 && $escape_counter < 100) {
	    $Flag = 1;
	    $RFlag = 0;
	    $escape_counter++;
	    $randfrag = int(rand($numfrags));
	    $start_forward = $fragend[$randfrag] - $fraglen[$randfrag] + 1;
	    $end = $fragend[$randfrag];
	    $len = $fraglen[$randfrag];
	    if($len >= $medfraglength && $len > ($medfraglength + $n * ($maxfraglength - $medfraglength))) {
		$Flag = 0;
	    }
	    if($len <= $medfraglength && $len < ($medfraglength - $n * ($medfraglength - $minfraglength))) {
		$Flag = 0;
	    }
	    if($len < $READLENGTH) {
		$RFlag = 1;
	    }
	}
	if($varcov =~ /\S/) {
	    undef @A;
	    $STRAND = $strand{$genes[$genenumber]};
	    if($STRAND ne $IVT_strand[$gene2ivt[$genenumber]]) {
		@A0 = split(/\t/,$IVT_data[$gene2ivt[$genenumber]]);
		my $N = @A0;
		for(my $i=0; $i<$N; $i++) {
		    $A[$i] = $A0[$N-($i+1)];
		}
	    } else {
		@A = split(/\t/,$IVT_data[$gene2ivt[$genenumber]]);
	    }
	    $dist_from_start = $start_forward / $seqlength;
	    $index = int($dist_from_start * @A);
	    if(rand() > $A[$index]) {
		$Flag = 0;
	    }
	    $gpd = $gene2percentdropoff[$gene2ivt[$genenumber]];
	    if(rand() > $gpd) {
		$Flag = 0;
	    }
	}

	if($RFlag == 0) {
	    $forward_read = substr($SEQ, $start_forward-1, $READLENGTH);
	    $start_reverse = $end - $READLENGTH + 1;
	    $reverse_read = substr($SEQ, $start_reverse-1, $READLENGTH);
	    if(length($forward_read) < $READLENGTH || length($reverse_read) < $READLENGTH) {
		$AFlag = 0;
	    }
	}
	else {
	    if(length($forward_read) < $len || length($reverse_read) < $len) {
		$AFlag = 0;
	    }
	    $forward_read = substr($SEQ, $start_forward-1, $len);
	    $start_reverse = $end - $len + 1;
	    $reverse_read = substr($SEQ, $start_reverse-1, $len);
	    $length2add = $READLENGTH - length($forward_read);
	    for(my $i=0; $i<$length2add; $i++) {
		my $base = &getrandombase();
		$forward_read = $forward_read . $base;
	    }
	    $length2add = $READLENGTH - length($reverse_read);
	    for(my $i=0; $i<$length2add; $i++) {
		my $base = &getrandombase();
		$reverse_read =  $base . $reverse_read;
	    }
	}
	if($Flag == 0) {
	    $AFlag = 0;
	}
    }

    $FLhash{$len}++;

    $rev_rev = reversecomplement($reverse_read);
    $fa = ">seq.$CNT" . "a\n$forward_read\n>seq.$CNT" . "b\n$rev_rev\n";

    # THE FOLLOWING GETS THE TRUE COORDINATES OF THE READ:

    #    ****   FORWARD READ   ****

    # fragment FORWARD read by its deletions and feed each piece separately to getcoords(),
    # then contcatenate coords together, also must correct the start and length of the FORWARD
    # read for the insertions (if any)

    if($RFlag == 0) {
	$RL = $READLENGTH;
    } else {
	$RL = $len;
    }
    $readlength = $RL;
    $checker_forward = 0;

    # this loop adjusts the start of the read for insertions/deletions upstream of the start
    $insertion_count = 0;
    undef @insertion_loc;
    undef @insertion_length;
    $indels_in_read_running_length=0;
    $insertion_at_beggining_of_read=0;
    $offset = 0;
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($INDELS[$ind][0] < $start_forward && $INDELS[$ind][0]+$INDELS[$ind][1]>=$start_forward && $INDELS[$ind][1] > 0) {
	    $insertion_loc[0] = 0;
	    $insertion_length[0] = $INDELS[$ind][0]+$INDELS[$ind][1]-$start_forward + 1;
	    $indels_in_read_running_length = $indels_in_read_running_length + $insertion_length[0];
	    $insertion_count++;
	    $insertion_at_beggining_of_read = $insertion_at_beggining_of_read + $insertion_length[0];
	    $readlength = $readlength - $insertion_length[0];
	}
	if(($INDELS[$ind][1] > 0 && ($INDELS[$ind][0] < $start_forward - 1)) || ($INDELS[$ind][1] < 0 && ($INDELS[$ind][0] <= $start_forward - 1))) {
	    if($INDELS[$ind][0] + $INDELS[$ind][1] < $start_forward) {
		$start_forward = $start_forward - $INDELS[$ind][1];
	    } else {
		$start_forward = $start_forward - $INDELS[$ind][1] + ($INDELS[$ind][0] + $INDELS[$ind][1] - $start_forward + 1);
	    }
	}
    }
    $end_forward = $start_forward + $readlength - 1;
    $coords = "";
    $cigar_string = "";
    # this loop does the fragmenting and feeding of each piece to getcoords()
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($start_forward <= $INDELS[$ind][0] && $INDELS[$ind][0] < ($start_forward+$readlength-1)) {
	    # get here if the indel is within the span of the read
	    if($INDELS[$ind][1] > 0) {  # insertion w.r.t. reference
                # in case insertion goes beyond end of read, don't want
                # to overcorrect later, the following 'if' takes care of that
		$ins_length = $INDELS[$ind][1];
		if( ($INDELS[$ind][0]+$INDELS[$ind][1]) > ($start_forward+$readlength-1) ) {
		    $adjustement_factor = (($INDELS[$ind][0]+$INDELS[$ind][1])-($start_forward+$readlength-1));
		    $readlength = $readlength + $adjustement_factor;  # yes, adding, so that subtraction later doesn't overcorrect
		    $end_forward = $end_forward + $adjustement_factor;
		    $checker_forward = 1;
		    $ins_length = $ins_length - $adjustment_factor;
		}
		if($start_forward <= $INDELS[$ind][0]) {
		    $insertion_loc[$insertion_count] = $INDELS[$ind][0] - $start_forward + 1 + $indels_in_read_running_length;
		    $insertion_length[$insertion_count] = $ins_length;
		    $indels_in_read_running_length = $indels_in_read_running_length + $ins_length;
		    $insertion_count++;
		}
	    }
	    if($INDELS[$ind][1] < 0) {  # deletion w.r.t. reference
		$deletion_length = -1 * $INDELS[$ind][1];
		$fragment_start = $start_forward + $offset;
		$fragment_length = $INDELS[$ind][0] - $fragment_start + 1;
		$SEQNAME = "seq." . $CNT . "a";
		$coords_to_add = getcoords($fragment_start, $fragment_length, $STARTS, $ENDS, $CHR, $SEQNAME);
		@CRD = split(/, / , $coords_to_add);
		for($crd=0; $crd<@CRD; $crd++) {
		    @SEG = split(/-/, $CRD[$crd]);
		    $len_seg = $SEG[1] - $SEG[0] + 1;
		    if($crd==0) {
			$cigar_string = $cigar_string . $len_seg . "M";
		    } else {
			$intronlen = $SEG[0] - $prev_end - 1;
			$cigar_string = $cigar_string . $intronlen . "N" . $len_seg . "M";
		    }	
		    $prev_end = $SEG[1];
		}
		$cigar_string = $cigar_string . $deletion_length . "D";
		$indels_in_read_running_length = $indels_in_read_running_length - $deletion_length;
		$coords = $coords . ", " . $coords_to_add;
		$offset = $offset + $fragment_length - $INDELS[$ind][1];
	    }
	    $readlength = $readlength - $INDELS[$ind][1];
	    $end_forward = $end_forward - $INDELS[$ind][1];
	}
    }
    $fragment_start = $start_forward + $offset;
    $fragment_length = $end_forward - $fragment_start + 1;
    $SEQNAME = "seq." . $CNT . "a";
    $coords_to_add = getcoords($fragment_start, $fragment_length, $STARTS, $ENDS, $CHR, $SEQNAME);
    @CRD = split(/, / , $coords_to_add);
    for($crd=0; $crd<@CRD; $crd++) {
	@SEG = split(/-/, $CRD[$crd]);
	$len_seg = $SEG[1] - $SEG[0] + 1;
	if($crd==0) {
	    $cigar_string = $cigar_string . $len_seg . "M";
	} else {
	    $intronlen = $SEG[0] - $prev_end - 1;
	    $cigar_string = $cigar_string . $intronlen . "N" . $len_seg . "M";
	}	
	$prev_end = $SEG[1];
    }
    $coords = $coords . ", " . $coords_to_add;
    $coords =~ s/^\s*,\s*//;
    $coords =~ s/\s*,\s*$//;
    $coords1 = $coords;
    $coords1 =~ s/.*\t//;

    # up to now the (forward) cigar string has no insertions, now going to put in all the insertions
    $cigar_string2 = "";
    $matchcount=0;  # this will keep a running total number of bases of the read accounted for by matches or insertions
    $cigar_string_temp = $cigar_string;
    $match_total_length = 0;    
    while($cigar_string_temp =~ /^(\d+)([^\d])/) {
	$num = $1;
	$type = $2;
	if($type eq 'M') {
	    $match_total_length = $match_total_length + $num;
	}
	$cigar_string_temp =~ s/^\d+[^\d]//;
    }
    while($cigar_string =~ /^(\d+)([^\d])/) {
	$num = $1;
	$type = $2;
	if($type eq 'M') {
	    $M1 = $num;
	    $current_length = 0;  # this keeps track of the length just in this match
	    for($ins=0; $ins<@insertion_loc; $ins++) {
		if($insertion_loc[$ins] >= $matchcount && $insertion_loc[$ins] <= $matchcount + $num) {
		    $M1 = $insertion_loc[$ins] - $current_length - $matchcount;
		    if($M1 > 0) {
			$cigar_string2 = $cigar_string2 . $M1 . "M";
		    }
		    if($insertion_loc[$ins] + $insertion_length[$ins] > $RL) {
			$insertion_length[$ins] = $insertion_length[$ins] - (($insertion_loc[$ins] + $insertion_length[$ins]) - $RL);
		    }
		    $cigar_string2 = $cigar_string2 . $insertion_length[$ins] . "I";
		    $current_length = $current_length + $M1;
		    $matchcount = $matchcount + $insertion_length[$ins];
		    $M1 = $num - $current_length;
		}
	    }
	    if($M1 > 0) {
		$cigar_string2 = $cigar_string2 . $M1 . "M";
		$current_length = $current_length + $M1;
	    }
	    $matchcount = $matchcount + $num;
	}
	if($type eq 'D' || $type eq 'N') {
	    $cigar_string2 = $cigar_string2 . $num . $type;
	}
	$cigar_string =~ s/^\d+[^\d]//;
    }

    $cigar_string_a = $cigar_string2;
    if($READLENGTH > $len) {
	my $L = $READLENGTH - $len;
	$cigar_string_a = $cigar_string_a . $L . "S";
    }

    #    ****   REVERSE READ   ****

    # fragment REVERSE read by its deletions and feed each piece separately to getcoords(),
    # then contcatenate coords together, also must correct the start and length of the REVERSE
    # read for the insertions (if any)


    if($RFlag == 0) {
	$RL = $READLENGTH;
    } else {
	$RL = $len;
    }
    $readlength = $RL;
    $checker_reverse = 0;

    # this loop adjusts the start of the read for insertions/deletions upstream of the start
    $insertion_count = 0;
    undef @insertion_loc;
    undef @insertion_length;
    $indels_in_read_running_length=0;
    $insertion_at_beggining_of_read=0;
    $offset = 0;
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($INDELS[$ind][0] < $start_reverse && $INDELS[$ind][0]+$INDELS[$ind][1]>=$start_reverse && $INDELS[$ind][1] > 0) {
	    $insertion_loc[0] = 0;
	    $insertion_length[0] = $INDELS[$ind][0]+$INDELS[$ind][1]-$start_reverse + 1;
	    $indels_in_read_running_length = $indels_in_read_running_length + $insertion_length[0];
	    $insertion_count++;
	    $insertion_at_beggining_of_read = $insertion_at_beggining_of_read + $insertion_length[0];
	    $readlength = $readlength - $insertion_length[0];
	}
	if(($INDELS[$ind][1] > 0 && ($INDELS[$ind][0] < $start_reverse - 1)) || ($INDELS[$ind][1] < 0 && ($INDELS[$ind][0] <= $start_reverse - 1))) {
	    if($INDELS[$ind][0] + $INDELS[$ind][1] < $start_reverse) {
		$start_reverse = $start_reverse - $INDELS[$ind][1];
	    } else {
		$start_reverse = $start_reverse - $INDELS[$ind][1] + ($INDELS[$ind][0] + $INDELS[$ind][1] - $start_reverse + 1);
	    }
	}
    }
    $end_reverse = $start_reverse + $readlength - 1;
    $coords = "";
    $cigar_string = "";
    # this loop does the fragmenting and feeding of each piece to getcoords()
    for($ind=0; $ind<@INDELS; $ind++) {  # indels are sorted by location
	if($start_reverse <= $INDELS[$ind][0] && $INDELS[$ind][0] < ($start_reverse+$readlength-1)) {
	    if($INDELS[$ind][1] > 0) {
                # in case insertion goes beyond end of read, don't want
                # to overcorrect later, the following 'if' takes care of that
		$ins_length = $INDELS[$ind][1];
		if( ($INDELS[$ind][0]+$INDELS[$ind][1]) > ($start_reverse+$readlength-1) ) {
		    $adjustement_factor = (($INDELS[$ind][0]+$INDELS[$ind][1])-($start_reverse+$readlength-1));
		    $readlength = $readlength + $adjustement_factor;  # yes, adding, so that subtraction later doesn't overcorrect
		    $end_reverse = $end_reverse + $adjustement_factor;
		    $checker_reverse = 1;
		    $ins_length = $ins_length - $adjustment_factor;
		}
		if($start_reverse <= $INDELS[$ind][0]) {
		    $insertion_loc[$insertion_count] = $INDELS[$ind][0] - $start_reverse + 1 + $indels_in_read_running_length;
		    $insertion_length[$insertion_count] = $ins_length;
		    $indels_in_read_running_length = $indels_in_read_running_length + $ins_length;
		    $insertion_count++;
		}
	    }
	    if($INDELS[$ind][1] < 0) {  # deletion w.r.t. reference
		$deletion_length = -1 * $INDELS[$ind][1];
		$fragment_start = $start_reverse + $offset;
		$fragment_length = $INDELS[$ind][0] - $fragment_start + 1;
		$SEQNAME = "seq." . $CNT . "b";
		$coords_to_add = getcoords($fragment_start, $fragment_length, $STARTS, $ENDS, $CHR, $SEQNAME);
		@CRD = split(/, / , $coords_to_add);
		for($crd=0; $crd<@CRD; $crd++) {
		    @SEG = split(/-/, $CRD[$crd]);
		    $len_seg = $SEG[1] - $SEG[0] + 1;
		    if($crd==0) {
			$cigar_string = $cigar_string . $len_seg . "M";
		    } else {
			$intronlen = $SEG[0] - $prev_end - 1;
			$cigar_string = $cigar_string . $intronlen . "N" . $len_seg . "M";
		    }	
		    $prev_end = $SEG[1];
		}
		$cigar_string = $cigar_string . $deletion_length . "D";
		$indels_in_read_running_length = $indels_in_read_running_length - $deletion_length;
		$coords = $coords . ", " . $coords_to_add;
		$offset = $offset + $fragment_length - $INDELS[$ind][1];
	    }
	    $readlength = $readlength - $INDELS[$ind][1];
	    $end_reverse = $end_reverse - $INDELS[$ind][1];
	}
    }
    
    $fragment_start = $start_reverse + $offset;
    $fragment_length = $end_reverse - $fragment_start + 1;
    $SEQNAME = "seq." . $CNT . "b";
    $coords_to_add = getcoords($fragment_start, $fragment_length, $STARTS, $ENDS, $CHR, $SEQNAME);
    @CRD = split(/, / , $coords_to_add);
    for($crd=0; $crd<@CRD; $crd++) {
	@SEG = split(/-/, $CRD[$crd]);
	$len_seg = $SEG[1] - $SEG[0] + 1;
	if($crd==0) {
	    $cigar_string = $cigar_string . $len_seg . "M";
	} else {
	    $intronlen = $SEG[0] - $prev_end - 1;
	    $cigar_string = $cigar_string . $intronlen . "N" . $len_seg . "M";
	}	
	$prev_end = $SEG[1];
    }
    $coords = $coords . ", " . $coords_to_add;
    $coords =~ s/^\s*,\s*//;
    $coords =~ s/\s*,\s*$//;
    $coords2 = $coords;
    $coords2 =~ s/.*\t//;
    $coords1 =~ /^(\d+)-/;
    $start_a = $1;
    $coords1 =~ /-(\d+)$/;
    $end_a = $1;
    $coords2 =~ /^(\d+)-/;
    $start_b = $1;
    $coords2 =~ /-(\d+)$/;
    $end_b = $1;
    $mergedcoords = merge($coords1, $coords2);
    if($mergedcoords eq "error") {
	$return_vector[0] = "none";
	$return_vector[1] = "none";
	$return_vector[2] = "none";
	return \@return_vector;
    }

    # up to now the (reverse) cigar string has no insertions, now going to put in all the insertions
    $cigar_string2 = "";
    $matchcount=0;  # this will keep a running total number of bases of the read accounted for by matches or insertions
    $cigar_string_temp = $cigar_string;
    $match_total_length = 0;    
    $cigar_string_temp = $cigar_string;
    $match_total_length = 0;    
    while($cigar_string_temp =~ /^(\d+)([^\d])/) {
	$num = $1;
	$type = $2;
	if($type eq 'M') {
	    $match_total_length = $match_total_length + $num;
	}
	$cigar_string_temp =~ s/^\d+[^\d]//;
    }
    while($cigar_string =~ /^(\d+)([^\d])/) {
	$num = $1;
	$type = $2;
	if($type eq 'M') {
	    $M1 = $num;
	    $current_length = 0;  # this keeps track of the length just in this match
	    for($ins=0; $ins<@insertion_loc; $ins++) {
		if($insertion_loc[$ins] >= $matchcount && $insertion_loc[$ins] <= $matchcount + $num) {
		    $M1 = $insertion_loc[$ins] - $current_length - $matchcount;
		    if($M1 > 0) {
			$cigar_string2 = $cigar_string2 . $M1 . "M";
		    }
		    if($insertion_loc[$ins] + $insertion_length[$ins] > $RL) {
			$insertion_length[$ins] = $insertion_length[$ins] - (($insertion_loc[$ins] + $insertion_length[$ins]) - $RL);
		    }
		    $cigar_string2 = $cigar_string2 . $insertion_length[$ins] . "I";
		    $current_length = $current_length + $M1;
		    $matchcount = $matchcount + $insertion_length[$ins];
		    $M1 = $num - $current_length;
		}
	    }
	    if($M1 > 0) {
		$cigar_string2 = $cigar_string2 . $M1 . "M";
		$current_length = $current_length + $M1;
	    }
	    $matchcount = $matchcount + $num;
	}
	if($type eq 'D' || $type eq 'N') {
	    $cigar_string2 = $cigar_string2 . $num . $type;
	}
	$cigar_string =~ s/^\d+[^\d]//;
    }
    $cigar_string_b = $cigar_string2;
    if($READLENGTH > $len) {
	my $L = $READLENGTH - $len;
	$cigar_string_b = $cigar_string_b . $L . "S";
    }

    $cigstr1 = "";
    $cigstr2 = "";
    if($rev_flip == 0) {
	$cigstr1 = $cigstr1 . "seq.$CNT";
	$cigstr1 = $cigstr1 . "a\t$CHR\t$start_a\t$end_a\t$len\t";
	$cigstr1 = $cigstr1 . "$cigar_string_a\t$coords1\t+\t";

	$cigstr2 = $cigstr2 . "seq.$CNT";
	$cigstr2 = $cigstr2 . "b\t$CHR\t$start_b\t$end_b\t$len\t";
	$cigstr2 = $cigstr2 . "$cigar_string_b\t$coords2\t-\t";
    } else {
	$cigstr1 = $cigstr1 . "seq.$CNT";
	$cigstr1 = $cigstr1 . "a\t$CHR\t$start_b\t$end_b\t$len\t";
	$cigstr1 = $cigstr1 . "$cigar_string_b\t$coords2\t-\t";

	$cigstr2 = $cigstr2 . "seq.$CNT";
	$cigstr2 = $cigstr2 . "b\t$CHR\t$start_a\t$end_a\t$len\t";
	$cigstr2 = $cigstr2 . "$cigar_string_a\t$coords1\t+\t";
    }

    @D = split(/, /,$mergedcoords);
    $bed = "";
    for($dd=0; $dd<@D; $dd++) {
	$D[$dd] =~ /(\d+)-(\d+)/;
	$SSS123 = $1;
	$EEE123 = $1;
	$D[$dd] =~ s/-/\t/;
	$str = "$chr{$GENE}\t$D[$dd]\t+\n";
	if($str =~ /^\S+\t\d+\t\d+\t\+$/ && $SSS123 <= $EEE123) {
	    if($seq_num_in_bedfile eq "true") {
		$bed = $bed . "seq.$CNT\t$str";
	    } else {
		$bed = $bed . "$str";
	    }
	}
	else {
	    print STDERR "-------\nsomething is wrong - check the log file\n";
	    print SIMLOGOUT "-------\nstr = $str (something is wrong)\n";
	    print SIMLOGOUT "readlength = $readlength\n";
	    print SIMLOGOUT "coords1 = $coords1\n";
	    print SIMLOGOUT "coords2 = $coords2\n";
	    print SIMLOGOUT "mergedcoords = $mergedcoords\n";
	    print SIMLOGOUT "3:GENE = $GENE\n";
	    print SIMLOGOUT "SEQ = $SEQ\n";
	    print SIMLOGOUT "$fa";
	    print SIMLOGOUT "STARTS = $STARTS\n";
	    print SIMLOGOUT "ENDS = $ENDS\n";
	    for($ind=0; $ind<@INDELS; $ind++) {
		print SIMLOGOUT "INDELS[$ind][0] = $INDELS[$ind][0]\n";
		print SIMLOGOUT "INDELS[$ind][1] = $INDELS[$ind][1]\n";
	    }
	    print SIMLOGOUT "start = $start\n";
	    print SIMLOGOUT "fragment end = $end\n";
	    print SIMLOGOUT "start_forward = $start_forward\n";
	    print SIMLOGOUT "end_forward = $end_forward\n";
	    print SIMLOGOUT "start_reverse = $start_reverse\n";
	    print SIMLOGOUT "end_reverse = $end_reverse\n";
	    print SIMLOGOUT "fragment_start = $fragment_start\n";
	    print SIMLOGOUT "fragment_length = $fragment_length\n";
	    print SIMLOGOUT "cutpoint = $cutpoint\n";
	    print SIMLOGOUT "seqlength = $seqlength\n";
	    print SIMLOGOUT "flip = $flip\n";
	    print SIMLOGOUT "checker_forward = $checker_forward\n";
	    print SIMLOGOUT "checker_reverse = $checker_reverse\n";
	    close(SIMLOGOUT);
	    open(SIMLOGOUT, ">>$logfilename");
	}
    }
    $return_vector[0] = $fa;
    $return_vector[1] = $bed;
    $return_vector[2] = $mergedspans;
    $return_vector[3] = $cigstr1;
    $return_vector[4] = $cigstr2;

    return \@return_vector;
}

sub merge () {
    ($aspans, $bspans) = @_;
    undef @astarts2;
    undef @aends2;
    undef @bstarts2;
    undef @bends2;

    $aspans =~ /(\d+)$/;
    $aspans_end = $1;
    $bspans =~ /^(\d+)/;
    $bspans_start = $1;
    if(!($aspans_end =~ /\S/) || !($bspans_start =~ /\S/)) {
	print SIMLOGOUT "Something is wrong:\naspans = $aspans\nbspans=$bspans\n\n";
	return "error";
    }

    @a = split(/, /, $aspans);
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$astarts2[$i] = $b[0];
	$aends2[$i] = $b[1];
    }
    @a = split(/, /, $bspans);
    for($i=0; $i<@a; $i++) {
	@b = split(/-/,$a[$i]);
	$bstarts2[$i] = $b[0];
	$bends2[$i] = $b[1];
    }
    if($aends2[@aends2-1] + 1 < $bstarts2[0]) {
	$merged_spans = $aspans . ", " . $bspans;
    }
    if($aends2[@aends2-1] + 1 == $bstarts2[0]) {
	$aspans =~ s/-\d+$//;
	$bspans =~ s/^\d+-//;
	$merged_spans = $aspans . "-" . $bspans;
    }
    if($aends2[@aends2-1] + 1 > $bstarts2[0]) {
	$merged_spans = $aspans;
	for($i=0; $i<@bstarts2; $i++) {
	    if($aends2[@aends2-1] >= $bstarts2[$i] && ($aends2[@aends2-1] <= $bstarts2[$i+1] || $i == @bstarts2-1)) {
		$merged_spans =~ s/-\d+$//;
		$merged_spans = $merged_spans . "-" . $bends2[$i];
		for($j=$i+1; $j<@bstarts2; $j++) {
		    $merged_spans = $merged_spans . ", $bstarts2[$j]-$bends2[$j]";
		}
	    }
	}
    }
    return $merged_spans;
}

sub reversecomplement () {
    ($sq) = @_;
    @A = split(//,$sq);
    $rev = "";
    for($i=@A-1; $i>=0; $i--) {
	$flag = 0;
	if($A[$i] eq 'A') {
	    $rev = $rev . "T";
	    $flag = 1;
	}
	if($A[$i] eq 'T') {
	    $rev = $rev . "A";
	    $flag = 1;
	}
	if($A[$i] eq 'C') {
	    $rev = $rev . "G";
	    $flag = 1;
	}
	if($A[$i] eq 'G') {
	    $rev = $rev . "C";
	    $flag = 1;
	}
	if($flag == 0) {
	    $rev = $rev . $A[$i];
	}
    }

    return $rev;
}

sub getcoords () {
    ($readstart, $readlength2, $starts, $ends, $CHR, $SEQNAME) = @_;
    $readend = $readstart + $readlength2 - 1;
    $COORDS = "";
    $starts =~ s/,\s*$//;
    $ends =~ s/,\s*$//;
    @S = split(/,/,$starts);
    @E = split(/,/,$ends);
    for($ii=0; $ii<@S; $ii++) {
	$S[$ii]++;
    }
    $cumlength = 0;
    $str = "";
    for($I=0; $I<@S; $I++) {
	$l = $E[$I] - $S[$I] + 1;
	$cumlength = $cumlength + $l;
	$str = $str . "$cumlength, ";
    }
    $prefix = 0;
    $exonnum = 0;
    while($prefix < $readstart) {
	$prefix = $prefix + ($E[$exonnum] - $S[$exonnum] + 1);
	$exonnum++;
    }
    $offset1 = $prefix - ($E[$exonnum-1] - $S[$exonnum-1] + 1);
    $firstexon = $exonnum - 1; # the read starts in this exon
    $prefix = 0;
    $exonnum = 0;
    while($prefix < $readend) {
	$prefix = $prefix + ($E[$exonnum] - $S[$exonnum] + 1);
	$exonnum++;
    }
    $offset2 = $prefix - ($E[$exonnum-1] - $S[$exonnum-1] + 1);
    $lastexon = $exonnum - 1; # the read ends in this exon
    if($lastexon > @S) {
	$lastexon = @S - 1;
    }
    for($en = $firstexon; $en<=$lastexon; $en++) {
	if($en == $firstexon) {
	    $Z = $S[$en] + ($readstart - $offset1 - 1);
	    $COORDS = "$Z";
	}
	else {
	    $Z = $S[$en];
	    $COORDS = $COORDS . ", $Z";
	    $W = $Z;
	    print SIMJUNCTIONSOUT "-$W\n";
	}
	if($en == $lastexon) {
	    $Z = $S[$en] + ($readend - $offset2 - 1);
	    $COORDS = $COORDS . "-$Z";
	}
	else {
	    $COORDS = $COORDS . "-$E[$en]";
	    $W = $E[$en];
	    $SEQNAME_temp = $SEQNAME;
	    if($rev_flip == 1) {
		if($SEQNAME_temp =~ /a/) {
		    $SEQNAME_temp =~ s/a/b/;
		} else {
		    $SEQNAME_temp =~ s/b/a/;
		}
	    }
	    print SIMJUNCTIONSOUT "$SEQNAME_temp\t$CHR:$W";
	}
    }
    return $COORDS;
}

sub getrandombase () {
    $x = int(rand(4));
    if($x == 0) {
	return "A";
    }
    if($x == 1) {
	return "C";
    }
    if($x == 2) {
	return "G";
    }
    if($x == 3) {
	return "T";
    }
}

# Subroutine to pad each intron with the rest of the exons in the gene
sub getpaddedintron () {
    ($geneid, $intron) = @_;
    $intron--;
    if($intron >= 0) {
	$INTRON = $gene2intron{$geneid}[$intron];
	$INTRON =~ /^[^:]+:(\d+)-(\d+)/;
	$intronstart = $1;
	$intronend = $2;
    }
    else {
	$INTRON = "X";
    }
    @aa = @{$gene2exon{$geneid}};
    $numexons = @aa;
    $SEQ = "";
    $indelcounter=0; # counts the number of indels per *gene*
    $offset = 0;
    for($jj=0; $jj<$numexons; $jj++) {
	if($strand{$geneid} eq "+") {
	    $ii = $jj;
	} else {
	    $ii = $numexons - $jj - 1;
	}
	$exon = $aa[$jj];
	@aa2 = @{$indels{$aa[$ii]}};
	# add to hash that associates the *gene that has the intron $INTRON* to indels in gene w/intron specific coords (start of gene = 1)
	if($intron == $ii-1 && $strand{$geneid} eq "+") {
	    $offset = $offset + $intronend - $intronstart + 1;
	}
	for($k=0;$k<@aa2;$k++) {
	    $geneWithIntron2indel_temp{$geneid}{$INTRON}[$indelcounter][0] = $offset + $aa2[$k][0];
	    $geneWithIntron2indel_temp{$geneid}{$INTRON}[$indelcounter][1] = $aa2[$k][1];  # the length of the indel
	    $TTT = $geneWithIntron2indel_temp{$geneid}{$INTRON}[$indelcounter][0];
	    $indelcounter++;
	}
	$aa[$ii] =~ /^[^:]+:(\d+)-(\d+)/;
	$exonstart = $1;
	$exonend = $2;
	$offset = $offset + $exonend - $exonstart + 1;
	if($intron == $ii-1 && $strand{$geneid} eq "-") {
	    $offset = $offset + $intronend - $intronstart + 1;
	}
	if($strand{$geneid} eq "-") {
	    $SEQ = $exonseq{$aa[$jj]} . $SEQ;
	}
	else {
	    $SEQ = $SEQ . $exonseq{$aa[$jj]};
	}
	if($intron == $jj) {
	    if($strand{$geneid} eq "-") {
		$SEQ = $intronseq{$INTRON} . $SEQ;
	    } else {
		$SEQ = $SEQ . $intronseq{$INTRON};
	    }
	}
    }
    $FLAG = 0;
    @aa = @{$geneWithIntron2indel_temp{$geneid}{$INTRON}};
    while($FLAG == 0) {
	$FLAG = 1;
	for($ii=0; $ii<@aa-1; $ii++) {
	    if($aa[$ii][0] > $aa[$ii+1][0]) {
		$temp = $aa[$ii][0];
		$aa[$ii][0] = $aa[$ii+1][0];
		$aa[$ii+1][0] = $temp;
		$temp = $aa[$ii][1];
		$aa[$ii][1] = $aa[$ii+1][1];
		$aa[$ii+1][1] = $temp;
		$FLAG = 0;
	    }
	}
    }
    for($ii=0; $ii<@aa; $ii++) {
	$geneWithIntron2indel{$geneid}{$INTRON}[$ii][0] = $aa[$ii][0];
	$geneWithIntron2indel{$geneid}{$INTRON}[$ii][1] = $aa[$ii][1];
	$TTT = $geneWithIntron2indel{$geneid}{$INTRON}[$ii][0];
    }
    undef %geneWithIntron2indel_temp;

    return $SEQ;
}

sub add_sequencing_error () {
    ($read, $error_rate) = @_;

    $Rlength = length($read);
    $num_errors = random_binomial(1, $Rlength, $error_rate);
    $Qstring = $read;
    $Qstring =~ s/./J/g;
    for($ie=0; $ie<$num_errors; $ie++) {
	$loc = int(rand($Rlength));
	$origbase = substr($read,$loc,1);
	if($origbase =~ /[AGCT]/) {
	    $B = getrandombase();
	    while($B eq $origbase) {
	        $B = getrandombase();
	    }
	    if(length($read) > $loc) {
		$Z = substr($read,$loc,1,$B);
		$r2 = int(rand(9.999));
		switch($r2) {
		    case 0  {$r2 = "!"}
                    case 1  {$r2 = "\""}
                    case 2  {$r2 = "#"}
                    case 3  {$r2 = "\$"}
                    case 4  {$r2 = "\%"}
                    case 5  {$r2 = "&"}
                    case 6  {$r2 = "'"}
                    case 7  {$r2 = "\("}
                    case 8  {$r2 = "\)"}
                    case 9  {$r2 = "*"}
	        }
	        $Z = substr($Qstring,$loc,1,$r2);
	    }
        }
    }
    
    $return_vector[0] = $read;
    $return_vector[1] = $Qstring;

    return \@return_vector;
}

sub add_error_to_tail () {
    ($read, $low_qual_tail_length, $percent_of_tails_that_are_low_qual, $quality_of_tail, $Qstring) = @_;

    @S4 = split(//,$read);
    for($ie=0; $ie<@S4; $ie++) {
       if($S4[$ie] eq "N") {
          $r2 = int(rand(9.999));
          switch($r2) {
  	      case 0  {$r2 = "!"}
	      case 1  {$r2 = "\""}
	      case 2  {$r2 = "#"}
	      case 3  {$r2 = "\$"}
	      case 4  {$r2 = "\%"}
	      case 5  {$r2 = "&"}
	      case 6  {$r2 = "'"}
	      case 7  {$r2 = "\("}
	      case 8  {$r2 = "\)"}
	      case 9  {$r2 = "*"}
	  }
	  $Z = substr($Qstring,$ie,1,$r2);
       }
     }

    $Rnum = rand(1);
    if($Rnum < $percent_of_tails_that_are_low_qual) {
	$Rlength = length($read);
	for($ie=$Rlength-$low_qual_tail_length; $ie<@S4; $ie++) {
	    $Rnum2 = rand(1);
	    if($Rnum2 > $quality_of_tail) { # the lower $quality_of_tail is, the more likely we are to substitute
		$origbase = substr($read,$ie,1);
		if($origbase =~ /[AGCT]/) {
		    $B = getrandombase();
    		    while($B eq $origbase) {
		        $B = getrandombase();
		    }
		    $Z = substr($read,$ie,1,$B);
		    $r2 = int(rand(9.999));
		    switch($r2) {
		        case 0  {$r2 = "!"}
    		        case 1  {$r2 = "\""}
		        case 2  {$r2 = "#"}
		        case 3  {$r2 = "\$"}
		        case 4  {$r2 = "\%"}
		        case 5  {$r2 = "&"}
		        case 6  {$r2 = "'"}
		        case 7  {$r2 = "\("}
		        case 8  {$r2 = "\)"}
		        case 9  {$r2 = "*"}
                    }
		}
		$Z = substr($Qstring,$ie,1,$r2);
	    }
	}
    }
    $return_vector[0] = $read;
    $return_vector[1] = $Qstring;

    return \@return_vector;
}

sub format_large_int () {
    ($int) = @_;
    @a = split(//,"$int");
    $j=0;
    $newint = "";
    $n = @a;
    for($i=$n-1;$i>=0;$i--) {
	$j++;
	$newint = $a[$i] . $newint;
	if($j % 3 == 0) {
	    $newint = "," . $newint;
	}
    }
    $newint =~ s/^,//;
    return $newint;
}
