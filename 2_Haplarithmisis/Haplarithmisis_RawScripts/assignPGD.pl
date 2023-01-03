#!/usr/bin/perl

use strict;
use warnings 'all';
use POSIX;
use Sys::Hostname;
use Getopt::Long;
use File::Basename;
use DBI;
use FindBin;

umask 007;

my $REQUEST_DIR = "/home/zamani/PGT/TEST/PGD_Requests";
my $PROCESSED_REQUEST_DIR = "/home/zamani/PGT/TEST/Processed_PGD_Requests";
my $SCRIPTS_DIR = $FindBin::Bin;


if (`find $REQUEST_DIR -maxdepth 1 -mindepth 1 -type f -name "PGD*"`) {
    my $tmpFile = `mktemp $TMP_DIR/test.XXXXXXXXX`; chomp($tmpFile);
    my @allFiles = `find $REQUEST_DIR -maxdepth 1 -mindepth 1 -type f -name "PGD*"`;

    for my $pgdFile(@allFiles) {
	chomp($pgdFile);

        if (-e  $pgdFile) { #check required because other assignTasks might have consumed it since "find" above
            system("mv $pgdFile $tmpFile");

            if (-e $tmpFile) { #fails if mv failed due to permission problems 
                my $filename = basename($pgdFile);
		my $OLD_REQUEST_FILENAME = "$PROCESSED_REQUEST_DIR/$filename";
                system("touch $OLD_REQUEST_FILENAME;  cp $tmpFile $OLD_REQUEST_FILENAME");
                system("\\rm $tmpFile");


    }

}
	

