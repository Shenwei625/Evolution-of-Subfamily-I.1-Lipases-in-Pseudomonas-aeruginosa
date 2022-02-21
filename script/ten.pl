#! usr/bin/perl
open LIST, '<dist2.txt' or die "Cannot write to ls.out:$!";
open STDOUT, '>format.txt' or die "Cannot write to ls.out:$!";
my @list = <LIST>;
foreach (@list) {
    s/\.[0-9]/ /m;
    s/_//;
    print;
}
