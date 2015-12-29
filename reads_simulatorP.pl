use FindBin qw($Bin);
use lib ("$Bin/lib", "$Bin/lib/perl5");
use YAML;

if(@ARGV < 4) {
    die "
Usage: reads_simulatorP.pl <num reads> <name> <num reads per chunk> <config file> <args>

I recommend 1,000,000 reads per chunk.  It takes about 3-4 hours for
a node to generate 1,000,000 read (pairs).

<args> are the usual args you would give to reads_simulator.pl

";
}

$num_reads = $ARGV[0];
$name = $ARGV[1];
$num_reads_per_chunk = $ARGV[2];
$config_file = $ARGV[3];

#my $r = {"submission_cmd" => "bsub", "submission_memory_flag" => "-M",
#    "memory_one" => 7000, "memory_two" => 8000};

#print YAML::Dump($r), "\n";
open my $fh, '<', $config_file or die "error opening $config_file: $!";
my $data = do { local $/; <$fh> };
#print $data, "\n";
my $config = Load($data);

#print STDERR "Lala $config->{submission_cmd} \n";
#die $config->{submission_cmd}, "\n";


$argstring = "";
$configstem = "temp";
for($i=3; $i<@ARGV; $i++) {
    $argstring = $argstring . " " . $ARGV[$i];
    if($ARGV[$i] eq "-configstem") {
	$configstem = $ARGV[$i+1];
    }
}
$type = "fa";
if($argstring =~ /outputfq/) {
    $type = "fq";
}
if($argstring =~ /-varcov +([^\s]+)/) {
    $varcov = $1;
}

if($num_reads % $num_reads_per_chunk == 0) {
    $numchunks = $num_reads / $num_reads_per_chunk;
    $size_lastchunk = $num_reads_per_chunk;
} else {
    $numchunks = int($num_reads / $num_reads_per_chunk) + 1;
    $size_lastchunk = $num_reads % $num_reads_per_chunk;
}

print STDERR "will run with $numchunks chunks\n";

$name_temp = $name . ".temp";

$path = `pwd`;
chomp($path);
$path2 = $path;

if($argstring =~ /-outdir +([^ ]+)/) {
    $path = $1;
}
if($path =~ /\/$/) {
    $path =~ s/.$//;
}
if($path2 =~ /\/$/) {
    $path2 =~ s/.$//;
}

$path =~ s!//!/!g;
$path2 =~ s!//!/!g;
print "path=$path\n";
print "path2=$path2\n";

open(OUT, ">$path/sim_temp.sh");

open(OUT2,">$path/$name_temp.txt");
print OUT2 "waiting\n";
close(OUT2);

print OUT "perl $path2/scripts/reads_simulator.pl 10 $name_temp $argstring\n";
print OUT "echo 'done' > $path/$name_temp.txt\n";
close(OUT);

`$config->{submission_cmd} $config->{submission_memory_flag} $config->{memory_two} bash $path/sim_temp.sh`;

print STDERR "writing initial configs\n";

$a = `cat $path/$name_temp.txt`;
until($a =~ /done/s) {
    sleep(10);
    $a = `cat $path/$name_temp.txt`;
}

open(OUT, ">sim_temp.sh");
$subsname = "simulated_reads_substitutions_" . $name . ".temp.txt";
$indelsname = "simulated_reads_indels_" . $name . ".temp.txt";
$altsname = "simulated_reads_transcripts_" . $name . ".temp.txt";
$ivtname = "IVT2GENE.txt";

# `cp /home/ggrant/pgfi_data/gpfs/ggrant/simulator/IVT2GENE.txt /home/ggrant/pgfi_data/gpfs/ggrant/simulator/vcrefseq7/`;

print STDERR "firing off chunks\n";

for($i=0; $i<$numchunks; $i++) {
    if($i == $numchunks - 1) {
	$numreads_thischunk = $size_lastchunk;
    } else {
	$numreads_thischunk = $num_reads_per_chunk;
    }
    $cntstart = $i * $num_reads_per_chunk + 1;
    $j = $i+1;
    $name_thischunk = $name . ".$j";
    $shell_thischunk = "$path/sim_temp" . ".$j" . ".sh";
    open(OUT, ">$shell_thischunk");
    $argstring =~ s/-configstem +[^ ]+//g;
    print OUT "perl $path2/scripts/reads_simulator.pl $numreads_thischunk $name_thischunk -configstem $configstem -usesubs $path/$subsname -useindels $path/$indelsname -usealts $path/$altsname -usepcd -usevarcov $path/$ivtname -sn -cntstart $cntstart $argstring\n";
    print OUT "echo 'done' >> $path/$name_temp.txt\n";
    close(OUT);
    `$config->{submission_cmd} $config->{submission_memory_flag} $config->{memory_one} bash $shell_thischunk`;
}

$a = `wc -l $path/$name_temp.txt`;
$a =~ /(^\d+)/;
$n = $1;
while($n < $numchunks+1) {
    $m = $n-1;
    if($m > 0) {
	if($m > 1) {
	    print STDERR "$m chunks done out of $numchunks\n";
	} else {
	    print STDERR "1 chunk done out of $numchunks\n";
	}
    }
    sleep(5);
    $a = `wc -l $path/$name_temp.txt`;
    $a =~ /(^\d+)/;
    $n = $1;
}
if($numchunks > 1) {
    print STDERR "$numchunks chunks done out of $numchunks\n";
} else {
    print STDERR "$numchunks chunk done out of $numchunks\n";
}

print STDERR "all chunks finished, going to merge them now.\n";

print STDERR "perl $path2/scripts/merge_simulated.pl $name $numchunks $path $type\n";
`perl $path2/scripts/merge_simulated.pl $name $numchunks $path $type`;

print STDERR "cleaning up...\n";

print STDERR "perl $path2/scripts/cleanup_simulator_chunkfiles.pl $name $numchunks $path $type\n";
`perl $path2/scripts/cleanup_simulator_chunkfiles.pl $name $numchunks $path $type`;

print STDERR "ok all done.\n";
