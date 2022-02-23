#! usr/bin/perl
open LIST, '<out_file' or die "Cannot read from out_file:$!";
open STDOUT, '>region.txt' or die "Cannot write to region.txt:$!";
chomp(my @list = <LIST>);
foreach (@list) {
    if (/(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        print "$2:$9-$10\n";
    }
}