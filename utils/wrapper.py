# Universal wrapper, by Ulf Ekstrom 2011
import os, subprocess, shutil, glob, signal, optparse, sys, shlex, tarfile, tempfile, zipfile

def main():
    usage = \
	"""%prog [options] 
Create a temporary directory, copy in files, execute one or several commands and copy back results.
"""
    parser =optparse.OptionParser(usage)

    parser.add_option('-v','--verbose', action='store_true', dest="verbose",
		     help='Print more information')
    parser.add_option('-t','--tmpdir', type='string', default=None,dest="tmpdir",
		     help='Base directory for temporary directory',metavar="PATH")
    parser.add_option('-T','--tmpfull', type='string', default=None,dest="tmpfull",
		     help='Full path to temporary directory',metavar="PATH")
    parser.add_option('-k','--keep-tmp', action='store_true', dest="keep_tmp",default=False,
		     help='Keep temporary directory')

    group = optparse.OptionGroup(parser, 'Input options')
    infiles=[]
    def infilecb(option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
	    if len(ff) != 2:
		print >> sys.stderr, "Wrapper: Illegal argument to",opt
		sys.exit(1)
	    infiles.append([ff[1],ff[0]])
	else:
	    infiles.append([value])
    def replacecb(option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
	    ff = [ff[0],'='.join(ff[1:])]
	    if len(infiles) == 0:
		print >> sys.stderr, "Wrapper:",option,"given before any -i/--infile was specified"
		sys.exit(1)
	else:
	    print >> sys.stderr, "Wrapper: Illegal argument to",opt," expected NAME=VALUE"
	    sys.exit(1)
	infiles[-1].extend(ff)
    group.add_option('-i','--infile', type='string', action='callback',callback=infilecb,
		     help='Copy file to TMPDIR/NAME',metavar='[NAME=]FILE')
    group.add_option('-r','--replace', type='string', action='callback',callback=replacecb,
		     help='Replace NAME by VALUE in the previous --infile.',metavar="NAME=VALUE")
    group.add_option('-x','--extract', type='string', action='append',default=[],dest="extract",
		     help='Extract archive FILE (tar, tgz, zip) in the temporary directory',metavar='FILE')
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Command options')
    commands=[]
    env = os.environ.copy()
    def definecb(option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
            env[ff[0]] = ff[1]
	else:
	    print >> sys.stderr, "Wrapper: Illegal argument to",opt," expected NAME=VALUE"
	    sys.exit(1)
    def commandcb(option, opt, value, parser):
	commands.append(shlex.split(value))
    def pendcb(option, opt, value, parser):
	if len(commands) == 0:
		print >> sys.stderr, "Wrapper:",option,"given before any -c/--command was specified"
		sys.exit(1)
	if opt in ['-a','--append']:
	    commands[-1].extend(shlex.split(value))
	else:
	    commands[-1] = shlex.split(value) + commands[-1]
    group.add_option('-c','--command', type='string', action='callback',callback=commandcb,
		     help='Add a new command to execute',metavar='COMMAND')
    group.add_option('--append', type='string', action='callback',callback=pendcb,
		     help='add ARGS last on the previous command line given by -c/--command',metavar='ARGS')
    group.add_option('--prepend', type='string', action='callback',callback=pendcb,
		     help='add ARGS first on the previous command line given by -c/--command',metavar='ARGS')
    group.add_option('-D','--define', type='string', action='callback',callback=definecb,
		     help='Define NAME=VALUE in the environment of the commands',metavar='NAME=VALUE')
    parser.add_option_group(group)

    outfiles=[]
    def getfilecb(option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
	    if len(ff) != 2:
		print >> sys.stderr, "Wrapper: Illegal argument to",opt
		sys.exit(1)
	    outfiles.append(ff)
	else:
	    outfiles.append([value])
    group = optparse.OptionGroup(parser, 'Output options')
    group.add_option('-g','--getfile', type='string', action='callback',callback=getfilecb,
		     help='Copy from to TMPDIR/NAME to FILE (NB argument order). Wildcards are allowed if no FILE is given.',
                     metavar='NAME[=FILE]')
    group.add_option('-o','--stdout', type='string', action='store',default=None,dest="stdout",
		     help='Write stdout to FILE',metavar="FILE")
    group.add_option('-e','--stderr', type='string', action='store',default=None,dest="stderr",
		     help='Write stderr to FILE',metavar="FILE")
    group.add_option('-a','--archive', type='string', action='store',default=None,dest="archive",
		     help='Path to output archive (tgz) file',metavar="FILE")
    group.add_option('-A','--archive-files', type='string', action='append',default=[],dest="archfiles",
		     help='Files matching this pattern will be achived',metavar="PATTERN")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(sys.argv[1:])

    try: #Create temporary directory
        if options.tmpfull is None:
            if options.tmpdir is None:
                tmppath = tempfile.mkdtemp()
            else:
                tmppath = tempfile.mkdtemp(dir=options.tmpdir)
            if options.verbose: print "Wrapper: Created temporary directory",tmppath
        else:
            tmppath = os.path.abspath(options.tmpfull)
            if os.path.isdir(tmppath):
                if options.verbose: print "Wrapper: Temporary directory '%s' already there" % tmppath
                options.keep_tmp = True
            else:
                if options.verbose: print "Wrapper: Creating Temporary directory '%s'" % tmppath
                os.makedirs(tmppath)
    except:
        print >> sys.stderr, "Wrapper: Cannot create temporary directory '%s', quitting." % tmppath
        sys.exit(1)
    try:
        if options.stdout is not None:
            try:
                stdout = open(options.stdout,"w")
            except:
                print >> sys.stderr, "Wrapper: Cannot create output file",options.stdout,"quitting."
                sys.exit(1)
        else:
            stdout = None
        if options.stderr is not None:
            try:
                stderr = open(options.stderr,"w")
            except:
                print >> sys.stderr, "Wrapper: Cannot create error file",options.stderr,"quitting."
                sys.exit(1)
        else:
            stderr = None
	run_wrapped(tmppath,commands,copyin=infiles,copyout=outfiles,verbose=options.verbose,env=env,
                    stdout=stdout,stderr=stderr,archive=options.archive,archfiles=options.archfiles,
                    extract=options.extract)
    except:
	print >> sys.stderr, "Wrapper: DEBUG: Something went wrong, re-raising Python error.."
	raise
    finally:  # Delete temporary directory
        if not options.keep_tmp:
            shutil.rmtree(tmppath)
            if options.verbose: print "Wrapper: Temporary directory deleted"
        else:
            if options.verbose: print "Wrapper: Temporary directory kept:",tmppath


def run_wrapped(tmppath, cmdlines, stdin=None, stdout=None, stderr=None,
                copyin=[], copyout=[], verbose=False, 
                env=None,extract=[],archive=None,archfiles=[]):
    """Create tmpdir, cd there and execute cmdlines, a list of arguments
    to Popen.
    cmdlines is a list of command lines in list form, ex: cmdlines=[['ls','-l'],[..]]
    copyin is a list of files to copy, in the form
    copyin=[(frompath,localname,rep1,with1,rep2,with2,..)] where rep1 will be 
    replaced by with1 during the copy.
    localname can be omitted, then the basename of frompath will be used.
    copyout is similar but no replacement is done. [(topath,localname),..]
    """
    for f in copyin: # Copy in files, stop if problem
        src = f[0]
        if len(f) % 2 == 0: # dst name given
            dst = os.path.join(tmppath,f[1])
            rep = f[2:]
        else:
            dst = os.path.join(tmppath,os.path.basename(src))
            rep = f[1:]
        rep = zip(rep[0::2],rep[1::2])
        try:
            copy_with_replacement(src,dst,rep,verbose=verbose)
        except:
            print >> sys.stderr, "Wrapper: File copy failed, quitting."
            raise
    # Extract archives, stop if problem
    def hassuffixes(a,ss):
        for s in ss:
            i = a.rfind(s)
            if i > 0 and i == len(a) - len(s):
                return True
        return False
    for a in extract:
        if verbose: print "Wrapper: Extracting archive",a
        if hassuffixes(a,['.tar']):
            t = tarfile.open(a,'r')
        elif hassuffixes(a,['.tgz','.tar.gz']):
            t = tarfile.open(a,'r:gz')
        elif hassuffixes(a,['.tbz','.tar.bz2']):
            t = tarfile.open(a,'r:bz2')
        elif hassuffixes(a,['.zip']):
            t = zipfile.ZipFile(a)
        else:
            print >> sys.stderr, "Wrapper: Cannot determine type of archive file",a,"quitting."
            sys.exit(1)
        try:
            t.extractall(path=tmppath)
        except:
            print >> sys.stderr, "Wrapper: Error extracting files from",a,"quitting."
            raise
    workdir = os.getcwd()
    try: # Run commands, stop if one fails
        os.chdir(tmppath)
        for cmd in cmdlines:
            pgid = None
            oldhand = []
            try:
                cmd = "'"+"' '".join(cmd)+"'"
                def handler(signum,frame):
                    if pgid is not None: 
                        if signum == signal.SIGALRM: signum = signal.SIGTERM
                        os.killpg(pgid,signum)
                catch = [signal.SIGTERM, signal.SIGHUP, signal.SIGINT, signal.SIGALRM]
                for s in catch:
                    oldhand.append(signal.signal(s,handler))
                if verbose:
                    print "Wrapper: Executing",cmd
                p = subprocess.Popen(cmd,stdin=stdin,stdout=stdout,stderr=stderr,shell=True,
                                     preexec_fn=os.setsid,env=env)
                pgid = p.pid
                p.wait() # OSError will be raised here if killed
                if p.returncode != 0:
                    raise subprocess.CalledProcessError(p.returncode,cmd)
            except subprocess.CalledProcessError, e:
                print >> sys.stderr, "Wrapper: Command "+cmd+" failed with returncode",e.returncode
                raise e
            except OSError, e: # Typically signal
                print >> sys.stderr, "Wrapper: Got OSError %i while executing" % e.errno,cmd
                raise e
            finally: # restore signal handlers
                if len(oldhand) > 0:
                    for s,h in zip(catch,oldhand):
                        signal.signal(s,h)
    finally: # Copy back files & create archive, always do this even if a command fails.
        # Create archive
        if archive is not None:
            try:
                archive = os.path.join(workdir,archive)
                if verbose: print "Wrapper: Creating archive",archive
                t = tarfile.open(archive,'w:gz')
                for p in archfiles:
                    for f in glob.glob(p):
                        t.add(f)
                t.close()
            except:
               print >> sys.stderr, "Wrapper: Error creating tar archive"
               raise
        # Copy back files
        os.chdir(workdir)
        def cpout(src,dst=None):
            if verbose:
                print "Wrapper: Copying back",src
            if os.path.isfile(src):
                if dst is None: dst='.'
                shutil.copy2(src,dst)
            elif os.path.isdir(src):
                if dst is None: dst=os.path.basename(src)
                if os.path.exists(dst): shutil.rmtree(dst)
                shutil.copytree(src,dst)
            elif verbose: # Just give a warning, don't panic
                print "Wrapper: When copying out, file",src,"not found."
        for f in copyout:
            if len(f) == 2:
                cpout(os.path.join(tmppath,f[0]),f[1])
            else: # allow patterns in copyout if there is no = 
                for ff in glob.glob(os.path.join(tmppath,f[0])):
                    cpout(ff)

def copy_with_replacement(src,dst,rep=[],verbose=False):
    try:
        if verbose: print "Wrapper: Copying '"+src+"' to '"+dst+"'"
        if len(rep) == 0:
            shutil.copy2(src,dst)
        else:
            if verbose:
                for r in rep:
                    print "Wrapper: Replacing '"+r[0]+"' by '"+r[1]+"'"
            of = open(dst,'w')
            for l in open(src,'r').readlines():
                for r in rep:
                    l = l.replace(r[0],r[1])
                of.write(l)
            of.close()
    except IOError,e:
        print >> sys.stderr, "Wrapper: Error copying file '"+src+"'"
        raise        

main()
