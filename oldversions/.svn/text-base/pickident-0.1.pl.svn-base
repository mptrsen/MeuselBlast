#!/usr/bin/perl
use strict;	# make me write clean code
use warnings;	# cry if something looks odd

use Bio::Seq;	# use sequence object stuff
use Bio::SeqIO;	# reading and writing sequence files

my $version = 0.1;

# handle command line args
use Getopt::Std;

my %opts;
getopts('1:2:3:4:l:t:q:', \%opts);
die "Usage: $0 -1 <inputfile1> -2 <inputfile2> -3 <inputfile3> [-l <logfile] [-n threshold] [-q query]\n" if (length(%opts) == 0);
my $inf_01	= $opts{1};
my $inf_02	= $opts{2};
my $inf_03	= $opts{3};
my $inf_04	= $opts{4};
my $logfile	= $opts{'l'};
my $threshold	= $opts{'t'};
my $query	= $opts{'q'};

my $search_string = 'Megajapyx';
my %read_ids_of_contig;
my %num_reads_of_contig;

# open logfile
my $logfh = &open_logfile($logfile);


# open file 01 and save contig IDs, length, #identities in separate arrays 
# (with indexes corresponding)
my @contig_ids;
my @lengths;
my @identities;

open(my $infh_01, '<', $inf_01) or die "wah! no open $inf_01 for reading: $!\n";

while (<$infh_01>) {
	chomp;
	if (/(\d+)\|$search_string/i) {
		push(@contig_ids, $1);
	}
	elsif (/Length = (\d+)/) {
		push(@lengths, $1);
	}
	elsif (/Identities = (\d+)/) {
		push(@identities, $1);
	}
}

close $infh_01;

# print out what we just collected, on the screen and into the logfile
print "collected " . @contig_ids . " contig IDs from $inf_01 (100% ident)\n";
print {$logfh} "collected " . @contig_ids . " contig IDs from $inf_01 (100% ident)\n";

# open file 02 and save in %read_ids_of_contig: contig ID => read IDs (as a single string)
# also save in %num_reads_of_contig: contig ID => #reads if #reads >= 10
open(my $infh_02, '<', $inf_02) or die "wah! no open $inf_02 for reading: $!\n";

# discard first line as it is a header we don't need
<$infh_02>;

my $count = 0;
my $found_count = 0;
while (<$infh_02>) {
	++$count;
	# pick out: contig ID, #reads, reads; save if #reads >=10
	/(\d+)\s+(\d+)\s+(.+)/;
	if ($2 >= $threshold) {
		$read_ids_of_contig{$1} = $3;
		$num_reads_of_contig{$1} = $2;
		++$found_count;
	}
}

close $infh_02;
printf {$logfh} "%d/%d contigs in %s where #reads >= %d\n", $found_count, $count, $inf_02, $threshold;
printf "%d/%d contigs in %s where #reads >= %d\n", $found_count, $count, $inf_02, $threshold;
print "\n";
sleep(1);

# print out the pretty table, push to @relevant_contig_ids those contigs that have:
# 	100% ident
#	#reads >= 10
my @relevant_contig_ids;
print "These contigs have identities of 100% and >=$threshold reads:\n\n";
print "ID\t\tlength\tidents\t#reads\tread IDs\n";
print '--------------------------------------------------------------------'."\n";
print {$logfh} "These contigs have identities of 100% and >=$threshold reads:\n\n";
print {$logfh} "ID\t\tlength\tidents\t#reads\tread IDs\n";
print {$logfh} '--------------------------------------------------------------------'."\n";
foreach my $i (1..$#contig_ids) {
	if (exists $read_ids_of_contig{$contig_ids[$i]}) {
		push(@relevant_contig_ids, $contig_ids[$i]);
		print $contig_ids[$i] . "\t" . $lengths[$i] . "\t" . $identities[$i] . "\t" . $num_reads_of_contig{$contig_ids[$i]} . "\t" . $read_ids_of_contig{$contig_ids[$i]} . "\n";
		printf {$logfh} "%d\t%d\t%d\t%d\t%s\n", $contig_ids[$i], $lengths[$i], $identities[$i], $num_reads_of_contig{$contig_ids[$i]}, $read_ids_of_contig{$contig_ids[$i]};
	}
}

print "\n";	# yay we did it
print {$logfh} "\n";

# now we're ready to skim through the mega file containing all the reads... 
# do this using SeqIO objects
printf "Going through %s seeking all reads for the relevant contigs...\n", $inf_03;
printf {$logfh} "Going through %s seeking all reads for the relevant contigs...\n", $inf_03;
my $seqio_obj = Bio::SeqIO->new(-file => $inf_03, -format => 'fasta');

$count = 1;
$found_count = 0;

# go through the fasta file one sequence at a time
while (my $seq_obj = $seqio_obj->next_seq) {
	# create a list of read IDs for each relevant contig
	foreach my $i (0..$#relevant_contig_ids) {
		my @list_of_reads = split(', ', $read_ids_of_contig{$relevant_contig_ids[$i]});
		# compare the read IDs with the headers (display_id) of the sequences in the fasta file
		foreach my $read_id (@list_of_reads) {
			# if the read ID is equal to the fasta header
			if ($read_id eq $seq_obj->display_id) {
				++$found_count;
				my $out_filename = '>>Contig'.$relevant_contig_ids[$i].'.fas';
				my $outf_obj = Bio::SeqIO->new(-file => $out_filename, -format => 'fasta');
				printf "read ID %s belongs to Contig%d (read #%d)\n", $read_id, $relevant_contig_ids[$i], $count;
				printf {$logfh} "read ID %s belongs to Contig%d (read #%d)\n", $read_id, $relevant_contig_ids[$i], $count;
				$outf_obj->write_seq($seq_obj);
			}
		}
	}
	++$count;
	if ($count % 1000 == 0) {
		printf "read ID %s... (%d reads so far)\n", $seq_obj->display_id, $count;
	}
}
print "found $found_count matches in $count sequences, collected into to the following files:\n";
print "Contig$_.fas\n" foreach (@relevant_contig_ids);
print {$logfh} "found $found_count matches in $count sequences, collected into to the following files:\n";
print {$logfh} "Contig$_.fas\n" foreach (@relevant_contig_ids);

sleep(3);
print "\n";
print {$logfh} "\n";

# okay we got the relevant reads in files with the contig IDs as names
# now find the contig in file 04 and append it to the correct file
printf "Going through %s picking out the relevant contigs...\n", $inf_04;
printf {$logfh} "Going through %s picking out the relevant contigs...\n", $inf_04;

$count = 1;
$found_count = 0;
my $contigs_file = Bio::SeqIO->new(-file => $inf_04, -format => 'fasta');

# go through the contigs file one seq at a time
while (my $contig = $contigs_file->next_seq) {
	foreach my $id (@relevant_contig_ids) {
		my $contig_id = $contig->display_id;
		$contig_id =~ /(\d+)/;
		if ($id eq $1) {
			++$found_count;
			my $out_filename = '>>Contig'.$id.'.fas';
			print "found a relevant contig ($id), appending to Contig$id.fas\n";
			print {$logfh} "found a relevant contig ($id), appending to Contig$id.fas\n";
			my $outf_obj = Bio::SeqIO->new(-file => $out_filename, -format => 'fasta');
			$outf_obj->write_seq($contig);
		}
	}
	printf "checking contig %s... (%d contigs so far)\n", $contig->display_id, $count if ($count % 1000 ==0);
	++$count;
}

printf "done. found %d contigs, appended to their respective files.\n", $found_count;
print "\n";
print "Thanks for using me!\n";
exit;


#####################
# functions follow
#####################

# open reference to logfile path
# input: filename
# output: file handle
sub open_logfile {
	my $logfile = shift @_;
	# open logfile
	open(my $logfh, '>', $logfile) or die "Wah: no open $logfile for writing: $!\n";
	printf {$logfh} "version $version\n";
	return $logfh;
}

# input: two-element list with contig ID and string of read IDs separated by ', '
# output: list of reads
sub get_list_of_reads_for_contig {
	my @contig_reads = @_;
	my @reads = split(', ' , $contig_reads[1]);
	return @reads;
}

# check if this is a valid FASTA file
# input: file name as string
# output: 1 if this is a valid FASTA file, 0 otherwise
sub validate_fasta_file {
	my $fastafile = shift @_;
	open(my $fh, $fastafile) or die "Wah: no open file for reading: $!\n";
	while (<$fh>) {
		next if /^\s+$/;
		return 0 if !/^[^>]/;
		return 1;
	}
	close $fh;
}
