use Data::Dumper;
use strict;
use Cwd;
# use Getopt::Long;
my $bowtie_index = '/u1/tools/public/Gmax_v9.0/Gmax_v_0.9';
my $gtf_file = '/u1/tools/public/Gmax_v9.0/annotation/Gmax_189_gene_cuff.gff3';
my $cuff_cmp = '/u1/tools/public/Gmax_v9.0/annotation/cuffcmp.combined.gtf';
my $fastqc = '/u1/tools/public/FastQC/fastqc';
if (scalar @ARGV != 1) {
	print "usage: perl rna_seq.pl config_file\n";
	exit;
} 
# GetOptions ("length=i" => \$length,    # numeric
              # "file=s"   => \$data,      # string
              # "verbose"  => \$verbose)   # flag
  # or die("Error in command line arguments\n");
chomp $ARGV[0];
my $config_file = $ARGV[0];
my %hash;
my @fastq = ();
read_config($config_file);
print Dumper(\%hash);
my @samples = grep/Sample.*rep.*/,sort keys %hash;
chdir($hash{'directory'});

my ($second, $minute, $hour) = localtime();
my $time="$hour:$minute:$second";

# open a writeable filehandle and print to the filehandle
open(my $outfile, '>', "$hash{'directory'}.log.txt");
print $outfile $time,"\n";


if($hash{'pipeline'} eq 'complete'){
foreach my $f(sort @samples){
push(@fastq,(grep/$hash{$f}/,glob("*.fastq")));
}
print $outfile join "\n",@fastq,"\n";
my $fastqc_com = $fastqc." -t ".scalar @fastq." ".join(" ",@fastq);
print $outfile $fastqc_com,"\n"; 
system "$fastqc_com 1>>$hash{'directory'}.log.txt 2>&1";
map {system "/u1/tools/public/trim_galore $_ 1>>$hash{'directory'}.log.txt 2>&1"} @fastq;
my $parent_dir = getcwd();
my %bam_files;
foreach my $read(@samples){
my ($trimfq) = grep/$hash{$read}/,glob("*.fq");
system "tophat  --bowtie1 -G $gtf_file  -o $hash{$read} -p 4 $bowtie_index $trimfq 1>>$hash{'directory'}.log.txt 2>&1";
chdir($hash{$read});
# system "mv accepted_hits.bam $hash{$read}.bam 1>>$hash{'directory'}.log.txt 2>&1";
print $outfile "printing flagstat for all reads bam file\n";
system "samtools flagstat accepted_hits.bam 1>>$hash{'directory'}.log.txt 2>&1";
my $header = "header_".$hash{$read};
system "samtools  view -H accepted_hits.bam >$header 1>>$hash{'directory'}.log.txt 2>&1";
#samtools view ./GSM915327_ILLUMINA/PAIRED/SRR486242/accepted_hits.bam|grep  NH:i:1 > uniq_GSM915327.bam"
my $bam = "$hash{$read}"."_unique.bam";
`samtools view accepted_hits.bam | grep NH:i:1 | cat $header - | samtools view -Sb - > $bam 1>>$hash{'directory'}.log.txt 2>&1`;
print $outfile "printing flagstat for unique reads bam file\n";
`samtools flagstat $bam 1>>$hash{'directory'}.log.txt 2>&1`; 
# system "rm $sam";
# system "mv $bam $hash{$read}.bam";
chdir($parent_dir);
my $cufflinks_read = $hash{$read}."/$hash{$read}"."_unique.bam";
$bam_files{$read} = $cufflinks_read;
my $cufflinks_outdir =  "cufflinks_".$hash{$read};
system "cufflinks -G $gtf_file -o $cufflinks_outdir -p 4 $cufflinks_read 1>>$hash{'directory'}.log.txt 2>&1";
system "rm $trimfq 1>>$hash{'directory'}.log.txt 2>&1";
}
my $output_dir = $hash{'Sample_1_name'}."_vs_".$hash{'Sample_2_name'};
system "cuffdiff --output-dir $output_dir $cuff_cmp -p 4 $bam_files{Sample_1_rep1},$bam_files{Sample_1_rep2} $bam_files{Sample_2_rep1},$bam_files{Sample_2_rep2} 1>>$hash{'directory'}.log.txt 2>&1";
}

elsif($hash{'pipeline'} eq 'cuffdiff_only'){
my %bam_files;
foreach my $read(@samples){
# print $read,"\n";
my $cufflinks_read = $hash{$read}."/$hash{$read}"."_unique.bam";
# print $cufflinks_read,"\n";
$bam_files{$read} = $cufflinks_read;
}
my $output_dir = $hash{'Sample_1_name'}."_vs_".$hash{'Sample2_name'};
print $outfile "cuffdiff --output-dir $output_dir $cuff_cmp -p 4 $bam_files{Sample_1_rep1},$bam_files{Sample_1_rep2} $bam_files{Sample_2_rep1},$bam_files{Sample_2_rep2}\n";
system "cuffdiff --output-dir $output_dir $cuff_cmp -p 4 $bam_files{Sample_1_rep1},$bam_files{Sample_1_rep2} $bam_files{Sample_2_rep1},$bam_files{Sample_2_rep2} 1>>$hash{'directory'}.log.txt 2>&1";
}

sub read_config{
my $file = shift;
my $counter = 0;
open(my $FH,$file);
while(<$FH>){
chomp;
next if(/\#/g);
my($key,$value) = split/\s+/,$_;
$hash{$key} = $value;
}
}
