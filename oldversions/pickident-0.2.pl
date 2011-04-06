#!/usr/bin/perl
=head1 NAME 

not defined yet, please baptise me :)

=head1 SYNOPSIS

pickident -s <search string> -B <BLAST result> -T <contig Trace list> [-R <Reads (raw data)>] [-C <Contigs (processed)>] [-t threshold] 

=head2 BLAST result file:

 the contig IDs for taxon <search string>
 sequence length
 # of identities and the percentage thereof (not yet implemented)

=head2 contig Trace file:

 number of reads
 list of read IDs

=head2 Reads FASTA file:

 all reads that correspond to each contig id if the #reads threshold is met or exceeded

=head2 Contigs FASTA file:

 the relevant contigs

... and then it will print out a pretty table in table.txt and one FASTA file for each relevant contig which will contain all the reads plus the contig at the end.

=head1 COPYRIGHT

(c) 2010 Malte Petersen <mptrsen@uni-bonn.de>

This code is provided as-is, i.e. without any warranty whatsoever. Use at your own risk. I take zero responsibility for anything you do with it, especially if something goes wrong.
=cut


use strict;	# make me write clean code
use warnings;	# cry if something looks odd

use Bio::Seq;	# use sequence object stuff
use Bio::SeqIO;	# reading and writing sequence files
use Benchmark;

my $t0 = new Benchmark;

my $version = 0.1;

# handle command line args
use Getopt::Std;

my $usage = '-s <search string> -B <BLAST result> -T <contig Trace list> [-R <Reads (raw data)>] [-C <Contigs (processed)>] [-t threshold]';

my %opts;
getopts('B:T:R:C:t:s:', \%opts);
die "Usage: $0 $usage\n" if (keys %opts == 0);
my $inf_01	= $opts{'B'};
my $inf_02	= $opts{'T'};
my $inf_03	= $opts{'R'};
my $inf_04	= $opts{'C'};
my $threshold	= $opts{'t'};
my $search_string = $opts{'s'};

die "Must define a search string (-s string)\n" if not defined $search_string;
chomp $search_string;
die "No BLAST result file provided\n" if not defined $inf_01;

my $logfile = 'log.txt';
my $tablefile = 'table.txt';
my @contig_ids;
my @lengths;
my @identities;
my %read_ids_of_contig;
my %num_reads_of_contig;
my %contig_refs;
my $count = 0;
my $found_count = 0;
my @relevant_contig_ids;
$threshold = 10 if not defined $threshold;

# print to a list of files simultaneously
# be careful not to clobber existing vars when using the 
# globtype syntax (*varname)!!
# also be careful with forking (|), danger of zombie invasion
#open( my $out, "|tee $logfile") or die "wah: $!\n";	# discarded this solution because it produces asynchroneous output
open(my $logfh, '>', $logfile) or die "Wah: no open $logfile for writing: $!\n";
print "Version $version\n";
print $logfh "Version $version\n";
print "Search string: $search_string\n";
print $logfh "Search string: $search_string\n";
print "Reads threshold: $threshold\n";
print $logfh "Reads threshold: $threshold\n";
print "\n";
print $logfh "\n";


# open file 01 and save contig IDs, length, #identities in separate arrays 
# (with indexes corresponding)

open(my $infh_01, '<', $inf_01) or die "wah! no open $inf_01 for reading: $!\n";

while (<$infh_01>) {
	++$count;
	chomp;
	if ( />[a-z]+(\d+)\|$search_string/i ) {
		push(@contig_ids, $1);
	}
	elsif (/Length = (\d+)/) {
		push(@lengths, $1);
	}
	elsif (/Identities = (\d+)\/(\d+)/) {
		push(@identities, $1);
	}
}

close $infh_01;


# print out what we just collected
print "collected " . @contig_ids . " contig IDs from $inf_01 (100% ident)\n\n";
print $logfh "collected " . @contig_ids . " contig IDs from $inf_01 (100% ident)\n\n";
if (scalar @contig_ids == 0) {
	print "no contig IDs found for search string '$search_string', exiting\n";
	print $logfh "no contig IDs found for search string '$search_string', exiting\n";
	exit;
}

# open file 02 and save in %read_ids_of_contig: contig ID => read IDs (as a single string)
# also save in %num_reads_of_contig: contig ID => #reads 
die "contig trace list not defined\n" if (!defined $inf_02);
open(my $infh_02, '<', $inf_02) or die "wah! no open $inf_02 for reading: $!\n";

# discard first line as it is a header we don't need
<$infh_02>;

while (<$infh_02>) {
++$count;
# pick out: contig ID, #reads, reads 
/(\d+)\s+(\d+)\s+(.+)/;
$num_reads_of_contig{$1} = $2;
$read_ids_of_contig{$1} = $3;
++$found_count;
}

close $infh_02;
printf "%d/%d contigs in %s where #reads >= %d\n", $found_count, $count, $inf_02, $threshold;
printf $logfh "%d/%d contigs in %s where #reads >= %d\n", $found_count, $count, $inf_02, $threshold;
print $logfh "\n";


# populate the references hash
foreach my $i (0..$#contig_ids) {
	my $id = $contig_ids[$i];
	$contig_refs{$contig_ids[$i]} = [$lengths[$i], $identities[$i]];
	# dereferencing the anonymous arrays in %contig_refs
	if (defined $num_reads_of_contig{$id}) {
		push( @{$contig_refs{$id}}, ( $num_reads_of_contig{$id}, $read_ids_of_contig{$id} ) );
	}
	else {
		push( @{$contig_refs{$id}}, ( 0, '') );
	}
}

# we don't need these lists anymore, so free memory
@contig_ids = ();
@lengths = ();
@identities = ();
%num_reads_of_contig = ();
%read_ids_of_contig = ();

# fuck this took me ages to figure out with all the reference 
# brainfuck
# sort values of %contig_refs by their 2nd element, 
# save corresponding keys in an array (ordered list)
my @sorted_keys = sort( { ${$contig_refs{$b}}[2] <=> ${$contig_refs{$a}}[2] } keys %contig_refs );


# print out the pretty table, push to @relevant_contig_ids those contigs that have:
# 	100% ident (granted)
#	#reads >= 10
&print_pretty_table;


print ( "tabular output printed into $tablefile\n" );	# phew we did it
print $logfh ( "\n" );
print $logfh ( "tabular output printed into $tablefile\n" );	# phew we did it
print $logfh ( "\n" );

my $t1 = new Benchmark;
my $td = timediff($t1, $t0);
print "time: ", timestr($td), "\n";

die "reads data file not defined\n" if (!defined $inf_03);

# now we're ready to skim through the mega file containing all the reads... 
# do this using SeqIO objects
printf "Going through %s seeking all reads for the relevant contigs...\n", $inf_03;
printf $logfh "Going through %s seeking all reads for the relevant contigs...\n", $inf_03;
my $seqio_obj = Bio::SeqIO->new(-file => $inf_03, -format => 'fasta');

$count = 1;
$found_count = 0;

# go through the fasta file one sequence at a time
while (my $seq_obj = $seqio_obj->next_seq) {
	# create a list of read IDs for each relevant contig
	foreach my $key (@relevant_contig_ids) {
		my @list_of_reads = split(', ', ${$contig_refs{$key}}[3]) if ${$contig_refs{$key}}[3];
		# compare the read IDs with the headers (display_id) of the sequences in the fasta file
		foreach my $read_id (@list_of_reads) {
			# if the read ID is equal to the fasta header
			if ($read_id eq $seq_obj->display_id) {
				++$found_count;
				my $out_filename = '>>Contig'.$key.'.fas';
				my $outf_obj = Bio::SeqIO->new(-file => $out_filename, -format => 'fasta');
				printf "read ID %s belongs to Contig%d (read #%d)\n", $read_id, $key, $count;
				printf $logfh "read ID %s belongs to Contig%d (read #%d)\n", $read_id, $key, $count;
				$outf_obj->write_seq($seq_obj);
			}
		}
	}
	++$count;
	if ($count % 1000 == 0) {
		printf "read ID %s... (%d reads so far)\n", $seq_obj->display_id, $count;
	}
}
print "\nfound $found_count matches in $count sequences, collected into to the following files:\n";
print $logfh "\nfound $found_count matches in $count sequences, collected into to the following files:\n";
foreach my $file (@relevant_contig_ids) {
	print "Contig$file.fas\n";
	print $logfh "Contig$file.fas\n";
}

sleep(3);
print "\n";
print $logfh "\n";

# okay we got the relevant reads in files with the contig IDs as names
# now find the contig in file 04 and append it to the correct file
die "Contig data file not defined\n" if (!defined $inf_04);
printf "Going through %s picking out the relevant contigs...\n", $inf_04;
printf $logfh "Going through %s picking out the relevant contigs...\n", $inf_04;

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
			print $logfh "found a relevant contig ($id), appending to Contig$id.fas\n";
			my $outf_obj = Bio::SeqIO->new(-file => $out_filename, -format => 'fasta');
			$outf_obj->write_seq($contig);
		}
	}
	last if ($found_count == scalar @relevant_contig_ids);
	printf "checking contig %s... (%d contigs so far)\n", $contig->display_id, $count if ($count % 1000 ==0);
	++$count;
}

printf "done. found %d contigs, appended to their respective files.\n", $found_count;
printf $logfh "done. found %d contigs, appended to their respective files.\n", $found_count;
print "\n";
print "Thanks for using me!\n";

exit;


#####################
# functions follow
#####################

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


# TODO: work on this!
sub print_pretty_table {

open( my $table, '>', $tablefile ) or die "wah: no open $tablefile: $!\n";

print $table "These contigs have identities of 100% and >=$threshold reads:\n\n";
print $table "ID\t\tlength\tidents\t#reads\tread IDs\n";
print $table '--------------------------------------------------------------------'."\n";
foreach my $key (@sorted_keys) {
	if ( exists $contig_refs{$key} and ${$contig_refs{$key}}[2] >= $threshold) {
		push( @relevant_contig_ids, $key );
		printf $table "%s\t%d\t%d\t%d\t%s\n",
			$key,	# contig ID (from @sorted_keys)
			${$contig_refs{$key}}[0],	# length
			${$contig_refs{$key}}[1],	# idents
			${$contig_refs{$key}}[2],	# #reads
			${$contig_refs{$key}}[3]	# read IDs (as string)
		;
	}
}
print $table "\n";
print $table "These contigs have identities of 100% and <$threshold reads:\n\n";
print $table "ID\t\tlength\tidents\t#reads\n";
print $table '--------------------------------------------------------------------'."\n";
foreach my $key (@sorted_keys) {
	if ( exists $contig_refs{$key} and ${$contig_refs{$key}}[2] < $threshold) {
		printf $table "%s\t%d\t%d\t%d\n",
			$key,	# contig ID (from @sorted_keys)
			${$contig_refs{$key}}[0],	# length
			${$contig_refs{$key}}[1],	# idents
			${$contig_refs{$key}}[2],	# #reads
		;
	}
}


close($table);

}
