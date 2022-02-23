#! usr/bin/perl
open LIST, '<dist2.txt' or die "Cannot read from dist2.txt:$!";
open STDOUT, '>format.txt' or die "Cannot write to format.txt:$!";
my @list = <LIST>;
foreach (@list) {
    s/\.[0-9]/ /m;
    s/_//;
    print;
}
