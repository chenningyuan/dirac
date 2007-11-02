#!/usr/bin/perl
#
#
# Usage: sfmakedepend [-s] [-e ext] [-f file] [-I incdir]
#                     [-m mod_ext] [-c] [-p] [-g] [-h] [-o obj_ext]
#                     [-a add_ext] [-pl print_level] [-d cpp_variables] file ...
#
#
# This is a makedepend script for Fortran, including Fortran 90.
# It searches for Fortran style includes, C preprocessor includes,
# and module dependencies to the extent that I understand them.
#
# Your files must have an extension listed in the @suffixes list
# below.  You might also want to modify $compile_string if you use
# the -c or -p option.  I call the compiler $(CFT) for historical
# reasons.
#
# The switch [-s] is for stupid Fortran compilers that don't know
# how to automatically send things through the C preprocessor.
# It is possible to force 'make' to invoke 'cpp' to create a .f
# file from a .F file (which has cpp directives), but make won't
# know that a .f file will depend on the files that the .F file
# included.  This option will provide those dependencies.
#
# The [-e ext] switch is used with the [-s] switch for compilers
# which expect an extension other than .f on source files.  For
# instance, for the Connection Machine one would use "-e fcm".
#
# The [-f file] switch is used to change the name of the current
# Makefile.
#
# The [-I incdir] option tells sfmakedepend to look in alternate
# directories for the include files.  There can be several "-I dir"
# options used at once.  The current directory is still searched
# first.
#
# The [-m mod_ext] option tells sfmakedepend what extension to
# use for Fortran 90 module files.  The default for "use My_Mod"
# is to list the dependency as "my_mod.mod" since this is what
# NAG f90 and IBM xlf both use.  Let me know if other compilers
#perl version is  use a different filename for the module information.
#
# The [-c] option specifies that the Cray compiler is being used.
# This compiler requires a "-p file.o" option for each object file
# that contains a module used by your current source file.
#
# The [-p] option specifies that the Parasoft compiler is being used.
# This compiler requires a "-module file.o" option for each object file
# that contains a module used by your current source file.
#
# The [-g] option specifies that the SGI compiler is being used.
# This compiler names the module file in uppercase with the extension
# .kmo.
#
# The [-h] option specifies that the HP compiler is being used.
# This compiler names the module file in uppercase with the extension
# .mod (added by Patrick Jessee who also fixed an include bug).
#
# The [-o obj_ext] option tells sfmakedepend what extension to use for
# object files.  The default is "o", but "obj", for instance, is
# appropriate on MS-DOS etc.  This option was added by Dave Love.
#
# The [-a add_ext] option (also added by Dave Love) tells sfmakedepend
# to add targets with extension add_ext to the rules for object files.
# For instance, to operate with (f77) ftnchek .prj files, you could use
# `-a prj' to get rules like:
# foo.prj foo.o: ...
#
# The [-i print_level] option contains the print level for printing
# debugging info of the script execution, default value is 0,what means
# no debugging print-out
#
#  The [-d cpp_variables] option is not implemented now.
#
#
# The final arguments contain the list of source files to be
# searched for dependencies.
#
# EXAMPLES
#       sfmakedepend -I /usr/local/include *.F
#       sfmakedepend -I -i 5 ../../include *.F90
#
# NOTES
#	This makedepend script is my first attempt at using perl 5
#	objects.  Therefore, it may not be the best example of how
#	to do this.  Also, it requires perl 5 and will die if you
#	to use it with an older perl.  The latest version is
#	available from:
#
#		http://marine.rutgers.edu/po/perl.html
#		ftp://ahab.rutgers.edu/pub/perl/sfmakedepend
#
#	Fortran 90 introduces some interesting dependencies.  Two
#	compilers I have access to (NAG f90 and IBM xlf) produce a
#	private "mod_name.mod" file if you define "module mod_name"
#	in your code.  This file is used by the compiler when you
#	use the module as a consistency check (type-safe).  On the
#	other hand, the Cray and Parasoft compilers store the module
#	information in the object file and then files which use the
#	modules need to be compiled with extra flags pointing to the
#	module object files.
#
#	This script assumes that all the files using and defining
#	modules are in the same directory and are all in the list of
#	files to be searched.  It seems that the industry has not
#	settled on a practical way to deal with a separate modules
#	directory, anyway.
#
#	I sometimes include non-existent files as a compile time
#	consistency check:
#
#	   #ifndef PLOTS
#	   #include "must_define_PLOTS"       /* bogus include */
#	   #endif
#
#	This program warns about include files it can't find, but
#	not if there is a "bogus" on the same line.
#
#     *	The f90 module dependencies can confuse some versions of
#	make, especially of the System V variety.  We use gnu
#	make because it has no problems with these dependencies.
#
# BUGS
#	It can sometimes produce duplicate dependencies.
#
#	It treats C preprocessor includes the same as Fortran
#	includes.  This can add unnecessary dependencies if you
#	use the -s flag and both kinds of includes.
#
#       Please let me know if you find any others.
#	Kate Hedstrom
#	kate@ahab.rutgers.edu
#	
######################################################################
#
#      recent changes by Miro Ilias, March 2007, Tel Aviv
#
# This script is still in the development in order to make it more 
# powerfull for exact Fotran90 dependencies.
#
######################################################################

package source_file;

# hashes containing names of included files, modules in files
%inc_files = ();
%main::mod_files = ();

# Constructor
sub new {
    my $type = shift;
    my $filename = shift;
    my $path = shift;
    my $self = {};
    $self->{'source_file'} = $filename;
    $self->{'filepath'} = $path;
    $self->{'includes'} = {};
    $self->{'uses'} = {};
    $self->{'modules'} = {};
    bless $self;
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub find_includes {
#==================================================================================
#  
#  
#            Search throughout the file to find keyword of included files
#
#   called recursively
#
#==================================================================================
    my $self = shift;
    my $file = $self->{'filepath'};
    my($after, $filepath, $ref, $included, $use, $modname);

    local(*FILE);
    local($_);

    if (-f $file) {
        open(FILE, $file) || warn "Can't open $file: $!\n";
    } elsif (-f "RCS/$file,v" || -f "$file,v" ) {
        system("co $file");
        open(FILE, $file) || warn "Can't open $file: $!\n";
        $main::rcs{$file} = 1;
    } else {
	return;
    }

   if ($main::iprint >= 5) {
     print "\tSEARCH throughout the open file: $file \n"; }

    while (<FILE>) {
	$included = "";
	$use = "";
        # look for Fortran style includes
        if (/^\s*include\s*['"]([^"']*)["']/i) {
            $included = $1; 
	    if ($main::iprint >=5) {
	    print "found \"include\"> $included \n"; }
	    $after = $';
	# C preprocessor style includes
        } elsif (/^#\s*include\s*["<]([^">]*)[">]/) {
            $included = $1;
	    if ($main::iprint >= 5) {
	      print "found #include> $included\n"; }
	    $after = $';
        # Fortran 90 "use"
        } elsif (/^\s*use\s+(\w+)/i) {   
	    $use = $1;
	    if ($main::iprint >=5) {
	    print "found use> $1\n"; }
# Make the module name lowercase except for SGI & HP.
# May be compiler dependent!!
	    if ($main::sgi || $main::hp) {
		$use = uc($use);
	    } else {
		$use = lc($use);
	    }
	    $self->{'uses'}{$use} = 1;
	### Fortran 90 module
	} elsif (/^\s*module\s+(\w+)/i) {
	    $modname = $1;
	#   $modname = "modules/".$1;
        #    print "MODNAME: $modname\n";
	    if ($main::iprint >=5 ) {
	      print "found module> $1\n"; }
	    if ($main::sgi || $main::hp) {
		$modname = uc($modname);
	    } else {
		$modname = lc($modname);
	    }
	    unless (lc($modname) eq "procedure") {
	        $main::mod_files{$modname} = $file;
		$self->{'modules'}{$modname} = 1;
	    }
	}

##############################################################################
####              ...  find file through "include"  ...                   ####
##############################################################################
	if ($included) {
	    if ( $inc_files{$included} ) {
		$filepath = $inc_files{$included}{'filepath'};
	    } else {
                $filepath = &main::findfile($included);
		$ref = new source_file($included, $filepath);
                $inc_files{$included} = $ref;
#              ....   Search included file for includes
		$ref->find_includes();
	    }
            if ( $filepath ) {
		$self->{'includes'}{$included} = 1;
	    } else {
                if ($after !~ /bogus/i) {
		    warn "Can't find file: $included\n";
		}
	    }
	}
    }
#   ...close file
    close FILE;
}
# end sub find_includes

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub print_includes {
#================================================================================
#
#
#    print include files (*.h, *.mod, *.inc, ...) into the makefile 
#  after print target
#
#
#  called recursively 
#
#================================================================================
    my $self = shift;
    my $target = shift;
    my $len_sum = shift;
    my($file, $ref);
    my %printed = ();

    if ($main::iprint >=5) {
        print STDOUT "called 'print_includes' for target >> $target:\n"; }

    foreach $file (keys %{$self->{'includes'}}) {

	next if $printed{$file};

	$ref = $inc_files{$file};

	my $len = length($ref->{'filepath'}) + 1;

	if (($len_sum + $len > 80) &&
	       (length($target) + 1 < $len_sum) ) {

	 #  if ($main::iprint >=5) { 
         #   print STDOUT "xxx len=$len ref=$ref ref-filepath=$ref->{'filepath'} :\n"; }

	    print "\n$target:";     # prints target into Makefile ...

            # control print out on the screen as well ...
	    if ($main::iprint >=5) { 
             print STDOUT "$target:\n"; }

	    $len_sum = length($target) + 1;

	  # if ($main::iprint >=5) { 
          #   print STDOUT "its len_sum:  $len_sum\n"; }

	}

     #  if ($main::iprint >=5) {
     #    print STDOUT "going to print...$ref->{'filepath'}  \n"; 
     #  }

	print " " . $ref->{'filepath'};

            # control print out on the screen
	    if ($main::iprint >=5) {
              if ( $ref->{'filepath'} )  {
               print STDOUT "* depends on " . $ref->{'filepath'} ."\n"; }
               else {print STDOUT "...empty string...\n"; }
            }

	$printed{$file} = 1;
	$len_sum += $len;

	$len_sum = $ref->print_includes($target, $len_sum);

	# if ($main::iprint >=5) { print STDOUT " " . $ref->print_includes($target, $len_sum);}
    }
    $len_sum;
}
# end sub print_includes

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub inc_mods {
#====================================================================
#
# return list of modules used by included files
#
#====================================================================
    my $self = shift;
    my($file, $ref, $mod, @sub_list);
    my @list = ();
    my %printed = ();

    foreach $mod (keys %{$self->{'uses'}}) {

        print STDOUT "mod is here ! $mod \n";

	print $mod,"\n";
	push(@list, $mod);
    }

    foreach $file (keys %{$self->{'includes'}}) {
	next if $printed{$file};
	$ref = $inc_files{$file};
	$printed{$file} = 1;
	@sub_list = $ref->inc_mods();
	@list = (@list, @sub_list);
    }
    @list;
#   print @list,"\n";
}
# end sub inc_mods

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub find_mods {
################################################################################
#
#  filenames containing the list of modules used by file and all its includes
#
################################################################################

    my $self = shift;
    my($ref, $modname, $file, @list, $base);
    my @module_files = ();
    my @mod_list = ();
    my @tmp_list = ();

# find modules used by include files
    if (%{$self->{'includes'}}) {
	foreach $file (keys %{$self->{'includes'}}) {
	    $ref = $inc_files{$file};
	    @list = $ref->inc_mods();
            @tmp_list = @mod_list;
	    @mod_list = (@tmp_list, @list);
	}
    }

# add them to the uses list (hash ensures uniqueness)
    foreach $modname (@mod_list) {
        $self->{'uses'}{$modname} = 1;
    }
    
#
# now find the filename that contains the module information
#

    foreach $modname (keys %{$self->{'uses'}}) {

#         $modname = "modules/".$modname;
         # print STDOUT "HERE modified : $modname\n";
         # print "HERE modified : $modname\n";
  #       print "HERE 1 ! xxxx/$modname\n";
  #       print " and $main::mymodnames\n";
  #
         if ( index($main::mymodnames,$modname)==-1) 
         { $modname="modules/".$modname; 
   #        print "module modname : $modname \n"
           }
         else
         {
             #  print " xxxx/$modname\n"; 
         }

	if ($main::cray || $main::parasoft) {
	    if ($file = $main::mod_files{$modname}) {
                $base = &main::basename($file, @main::suffixes);
		$file = $base . "." . $main::obj_ext;
		push(@module_files, $file);
	    } else {
		warn "Don't know where module $modname lives.\n";
	    }
	} else {
	    $modname .= "." . $main::mod_ext;
	    push(@module_files, $modname);
	}

  #       print "HERE 2 ! xxxx/$modname\n";
  #       print " and $main::mymodnames\n";
  #       if ( index($main::mymodnames,$modname)==-1) 
  #       { $modname="modules/".$modname; 
  #         print "module modname : $modname \n"
  #         }
  #       else
  #       { print " xxxx/$modname\n"; 
  #       }

    }
    sort(@module_files);
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub print {
#===========================================================================
#
#
#
#         print basic (include,modules) dependencies of object files
#
#
#
#===========================================================================
    my $self = shift;
    my $source = $self->{'source_file'};
    my $compile_string = "\t" . '$(CFT) $(FFLAGS) -c';
    my ($base, $object, $modname, $flag, $target, $ftarget);

    $base = &main::basename($source, @main::suffixes);

    ###  ... define target
    $target = $base . "." . $main::obj_ext;

    if ($main::iprint >= 5) {
      print STDOUT "source=$source base=$base target=$target main::suffixes=@main::suffixes";
    }

    if ($main::stupid) {
	$ftarget = $base . "." . $main::ext;
    }

    if ($main::cray) {
	$flag = " -p ";
    } elsif ($main::parasoft) {
	$flag = " -module ";
    }

#-----------------------------------------------------
#         print out "include" dependencies
#-----------------------------------------------------
    if (%{$self->{'includes'}}) {
	my $len_sum = length($target) + 1;
	if ($main::add_ext) {
	    $target .= " $base.$main::add_ext";
	    $len_sum += length($base) + length($main::add_ext) + 2;
	}

        # ... print here ...
        print "$target:";
         $self->print_includes($target, $len_sum);
        print "\n";

       #if ($main::iprint >=5) {
       #  print STDOUT "...again target, $target: len_sum=$len_sum self is $self->{'source_file'} \
       #      $self->{'filepath'}  $self->{'includes'} $self->{'uses'}  $self->{'modules'} \n";
       #}

        if ($main::stupid) {
	    $len_sum = length($ftarget) + 1;
            print "$ftarget:";
            $self->print_includes($ftarget, $len_sum);
            print "\n";
        }
    }

# ---------------------------------------------------------
#    clean out "use" of modules in own file
# ---------------------------------------------------------
    foreach $mod ( keys %{$self->{'uses'}} ) {
	if ( ${$self->{'modules'}}{$mod} ) {
	    delete ${$self->{'uses'}}{$mod};
	}
    }
# ----------------------------------------------------------
#              print out "use" dependencies
# ----------------------------------------------------------
    if (%{$self->{'uses'}} || %{$self->{'includes'}}) {
	@module_files = $self->find_mods();
	my $len_sum = length($target) + 1;
     ##  print "volaco tu je 1: $target\n";
        if ($main::iprint >=5) {
         print STDOUT " ... again print out target: $target for modules ... \n";
         if (@module_files) {print STDOUT "...there are associated  modules \n";}
         else {print STDOUT "...no associated  modules !\n";}
        }
        # ... print target only if associated with modules !
        if (@module_files) {  
	 print "$target:"; }

	foreach $file (@module_files) {
	    my $len = length($file) + 1;
	    if (($len_sum + $len > 80) &&
	     (length($target) + 1 < $len_sum)) {
		print "\n$target:";
		$len_sum = length($target) + 1;
	    }
	    $len_sum += $len;
	    print " " . $file;
	    if ($main::iprint >=5) {
	      print STDOUT "module file: $file ... its source file (hash): $main::mod_obj{$file} \n";
            }
	    # ... adding object dependency through hash !
            print " " . $main::mod_obj{$file};
	}
	if ($main::need_f) {
	    my $len = length($ftarget) + 1;
	    if (($len_sum + $len > 80) &&
		   (length($target) + 1 < $len_sum)) {
		print "\n$target:";
		$len_sum = length($target) + 1;
	    }
	    print " " . $ftarget;
	}
	print "\n";
# extra Cray / Parasoft stuff
        if ($main::cray || $main::parasoft) {
	    print $compile_string;
	    foreach $file (@module_files) {
		print $flag . $file;
	    }
	    if ($main::stupid) {
	        print " " . $ftarget . "\n";
	    } else {
	        print " " . $source . "\n";
	    }
	}
    }
}
#end sub print
#end package source_file
#
#######################################
###    Start of the main program    ###
#######################################
package main;

check_perl_ver();

get_inline_args();

get_mystring($mystring);

# Search for the includes in all the files
search_includes();

# Create new Makefile with new dependencies.
open(MFILE, $mf) || die "can't read Makefile $mf: $!\n";
open(NMFILE, "> Makefile.new") || die "can't write Makefile.new: $!\n";
select(NMFILE);

while (<MFILE>) {
    if (!/$mystring/) {
      # print "pred subr_print\n";
        print;
    } else {
        last;
    }
}

 print $mystring, "\n";

##################################################################
#                                                                #
#               print out module dependencies                    #
#                                                                #
##################################################################
  if ($main::iprint >=5) {
     print STDOUT "\n======  print of module dependencies  ==========\n";
  }

#  first get list of module files (*.mod) arising from entering source (*.F90) files
$mymodnames="";
foreach $modname (sort keys(%mod_files)) {
 $mn=$modname.".".$mod_ext;
 $mymodnames=$mymodnames.$mn." ";
}

if ($main::iprint >=5) {
    print STDOUT "   locale module files: $mymodnames\n";
}


%main::mod_obj = ();

if ( !( $cray || $parasoft) ) {
   foreach $modname (sort keys(%mod_files)) {
        ($name, $path, $suffix) =
          &fileparse($sources{$mod_files{$modname}}->{'filepath'}, @suffixes);

	  $object = $path . "/" . $name . "." . $obj_ext;

	  $object =~ s#^\.//##; # convert .//xxx.o to xxx.o

	  $modnameindx = $modname.".".$mod_ext;
          $main::mod_obj{$modnameindx} = $object;

	if ($main::iprint >=5) {
	  print STDOUT "$modname.$mod_ext: $object \
	   ... in hash table: $main::mod_obj{$modnameindx} \n";
	 }
	  print "$modname.$mod_ext: $object\n";
    }
}

########################################################################
#
#
#       Now print out include and use dependencies in sorted order.
#
#
########################################################################
  if ($main::iprint >=5) {
     print STDOUT "\n======  print of object dependencies  ==========\n";
  }

foreach $target (sort keys(%sources)) {
  if ($main::iprint >=5) {
    print STDOUT "---> calling print for the target (source file) : $target \n";
  }
  $sources{$target}->print();
}

  sort_make();

#########################################################
#    							#
#             End of main function                      #
#    							#
#########################################################

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub sort_make {
# ----------------------------------------------------------------
#                    Sort out the Makefiles
# ----------------------------------------------------------------
rename($mf, "$mf.old") || warn "can't overwrite $mf.old: $!\n";
rename('Makefile.new', $mf) ||
     warn "can't move Makefile.new to $mf: $!\n";

# Delete those RCS files we checked out
foreach $file (keys %rcs) {
    unlink($file); }

} #end sub sort_make

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub get_mystring {
# ----------------------------------------------------------------
#
# ----------------------------------------------------------------
#            place in Makefile
#$mystring = '# DO NOT DELETE THIS LINE - used by make depend';
 $mystring = '# DO NOT DELETE THIS LINE -- make depend depends on it.';

}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub search_includes{
# ----------------------------------------------------------------
# Search for the includes in all the files
# ----------------------------------------------------------------
foreach $file (@ARGV) {
    if ($iprint >= 5)
     { print "\n\t >>>> source file in the argument list: $file \n"; }
    $sources{$file} = new source_file($file, $file);
 #   ... process the file
    $sources{$file}->find_includes();
 }
} #end sub search_includes

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub check_perl_ver{
#============================================================================
#
# check the current perl version 
#
#============================================================================
if ($] < 5.000) { die "Exit ! Need perl 5.000 or newer !!!!\n";}
else { print STDOUT "current perl version $] is OK...\n";}
}
#end sub check_perl_ver

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub get_inline_args{
#============================================
#
#        get inline arguments
#
#============================================
use File::Basename;
use Getopt::Long;

# list of supported suffixes
@suffixes = qw( .c .C .cc .cxx .cpp .f .F .fcm .FCM .f90 .F90 .for .head .struct);
# print "supported ".@suffixes." suffixes: @suffixes \n";

###   http://www.perl.com/doc/manual/html/lib/Getopt/Long.html 
GetOptions("s", "e=s", "f=s", "I=s@", "d=s", "m=s", "c", "p", "pl=i", "g", "h", "o=s", "a=s")
  || die "problem in GetOptions !" ;

# For compilers that don't invoke cpp for you
if ($opt_s) {
	$stupid = 1;
}
if ($opt_e) {
	$ext = $opt_e;
} else {
	$ext = "f";
}

# list of directories to be searched for include files, starting with current directory
if (@opt_I) {
    @incdirs = @opt_I;
} elsif (@opt_i) {
    @incdirs = @opt_i;
}
#print "include directories: @incdirs \n";

# list of -Dxxx preprocessor variables for which *.F, *.F90 source files are going to be checked
if ($opt_d) {
  $d_list = $opt_d;
#  print "d_list: $d_list \n";
   preproc_list($d_list) ;
#  print "NEW d_list: @d_list \n";
}

if ($opt_f) {
    $mf = $opt_f;
} elsif (-f "makefile") {
    $mf = 'makefile';
} else {
    $mf = 'Makefile';
}
# print "makefile is set to '$mf' \n";

if ( !(-f $mf)) {
    system "touch $mf";
}

#    extension used for compiler's private module information
if ($opt_m) {
  $mod_ext = $opt_m;
} else {
  $mod_ext = 'mod';
}

if ($opt_c) {
    $cray = 1;
}

if ($opt_p) {
    $parasoft = 1;
}

if ($opt_g) {
    $sgi = 1;
    $mod_ext = 'kmo';
}

if ($opt_h) {
    $hp = 1;
}

if ($opt_pl) {
    $iprint = $opt_pl; # assign the print level 
#    print "setting print level $opt_pl \n";
} else { $iprint = 0; }

  print STDOUT "debugging print level set to: $iprint \n";

#
# need to add some more dependencies so the .f file gets created
#
if ($stupid && ($cray || $parasoft)) {
    $need_f = 1;
}

if ($opt_c && $opt_p) {
 die "Doesn't make sense to have both Cray and Parasoft options!";
}

#  ... object file extension
if ($opt_o) {
    $obj_ext = $opt_o;
} else {
    $obj_ext = 'o';
}

# extension for additional targets (like .prj)
if ($opt_a) {
    $add_ext = $opt_a;
}

# 
# --- now print out entered arguments
# 
if ($iprint >= 5){
 print "\n *** Verbose output from this perl script ***\n";
 print "print level: $iprint\n";
 if (@opt_I || @opt_i) {
   print "read include directories : @incdirs \n"; }
 if ($opt_d) {
   print "preprocessor statements: $d_list\n"; }
 print "the makefile to be modified: $mf\n";
 print "Fortran 90 modules extension: $mod_ext\n";
 print "object files extension: $obj_ext\n";
 if ($opt_a) {
   print "additional target extension: $add_ext \n";
 }
}
}
# end sub ....

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub preproc_list {
#================================================================================
#  process list of preprocessor statements which are in the array @d_list
#================================================================================
  print "  ...  d_list:  $d_list \n";
  @d_list = split(/ /,$d_list);
  print "  splitted  @d_list \n";
  $i=0;
  foreach $D_var (@d_list) {
  # print "D_var in array: $D_var $i @d_list[$i] \n";
    $D_var = substr $D_var, 2;
    @d_list[$i]=$D_var;
   print "new member: @d_list[$i]\n ";
      $i++;
  }
  # print "newly modified array:  @d_list \n";
}

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
sub findfile {
#=============================================================================
#
# Let's see if we can find the included file.  Look in current
# directory first, then in directories from -I arguments.  Finally,
# look for RCS files in current directory.
#
#=============================================================================
    my $file = shift;
    my($found, $i, $filepath);
    $found = 0;
    if ( -f $file ) {
        $found = 1;
        $file =~ s#^\./##; # convert ./foo.h to foo.h
        return $file;
    }
    foreach $i (0 .. $#incdirs) {
        $filepath = $incdirs[$i]."/".$file;
        if ( -f $filepath ) {
	    $found = 1;
            $filepath =~ s#^\./##;  # convert ./foo.h to foo.h
            return $filepath;
        }
    }
#
#    see if it is a checked-in RCS file
#
    if (-f "RCS/$file,v" || -f "$file,v" ) {
        $found = 1;
        system("co $file");
	$filepath = $file;
        $rcs{$file} = 1;
    }
    if ( ! $found ) {
	$filepath = "";
    }
    $filepath;
}
#end sub findfile

#**********************************************************************
#                   End of the program file
#**********************************************************************
