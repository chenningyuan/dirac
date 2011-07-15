# Universal wrapper, by Ulf Ekstrom 2011
import os, subprocess, shutil, glob, signal, optparse, sys, shlex, tarfile, tempfile, zipfile

def wrap_dirac():
    usage = """%prog [options] 
Creates a temporary directory, copy in files, execute one or several 
commands and copy back results.
"""
    W = Wrapper()
    parser = optparse.OptionParser(usage)

    parser.add_option('-v','--verbose', action='store_true', dest="verbose",
		     help='Print more information')

    group = optparse.OptionGroup(parser, 'Dirac specific options')
    group.add_option('--mb', type='string', default=None, dest="mb",
                     help='Amount of memory for Dirac, in MB',metavar="MEGABYTES")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Temporary directory options')
    group.add_option('-t','--tmpdir', type='string', default=None,dest="tmpbase",
		     help='Base directory for temporary directory',metavar="PATH")
    group.add_option('-T','--tmpfull', type='string', default=None,dest="tmppath",
		     help='Full path to temporary directory. Will be created if not there.',metavar="PATH")
    group.add_option('-k','--keep-tmp', action='store_true', dest="keeptmp",default=False,
		     help="Keep temporary directory, default if it's preexisting.")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Input options')
    group.add_option('-i','--infile', type='string', action='callback',callback=W.infile_callback,
		     help='Copy file to TMPDIR/NAME',metavar='[NAME=]FILE')
    group.add_option('-r','--replace', type='string', action='callback',callback=W.replace_callback,
		     help='Replace NAME by VALUE in the previous --infile.',metavar="NAME=VALUE")
    group.add_option('-x','--extract', type='string', action='append',default=[],dest="extract",
		     help='Extract archive FILE (tar, tgz, zip) in the temporary directory',metavar='FILE')
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Command options')
    group.add_option('-c','--command', type='string', action='callback',callback=W.command_callback,
		     help='Add a new command to execute',metavar='COMMAND')
    group.add_option('--append', type='string', action='callback',callback=W.modcmd_callback,
		     help='add ARGS last on the previous command line given by -c/--command',metavar='ARGS')
    group.add_option('--prepend', type='string', action='callback',callback=W.modcmd_callback,
		     help='add ARGS first on the previous command line given by -c/--command',metavar='ARGS')
    group.add_option('-D','--define', type='string', action='callback',callback=W.define_callback,
		     help='Define NAME=VALUE in the environment of the commands',metavar='NAME=VALUE')
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, 'Output options')
    group.add_option('-g','--getfile', type='string', action='callback',callback=W.getfile_callback,
		     help='Copy from to TMPDIR/NAME to FILE (NB argument order). Wildcards are allowed if no FILE is given.',
                     metavar='NAME[=FILE]')
    group.add_option('-o','--stdout', type='string', action='store',default=None,dest="stdout",
		     help='Write stdout to FILE',metavar="FILE")
    group.add_option('-e','--stderr', type='string', action='store',default=None,dest="stderr",
		     help='Write stderr to FILE',metavar="FILE")
    group.add_option('--stdin', type='string', action='store',default=None,dest="stdin",
		     help='Read stdin from FILE',metavar="FILE")
    group.add_option('-a','--archive', type='string', action='store',default=None,dest="archive",
		     help='Path to output archive (tgz) file',metavar="FILE")
    group.add_option('-A','--archive-files', type='string', action='append',default=[],dest="archpatterns",
		     help='Files matching these pattern(s) will be achived',metavar="PATTERN")
    parser.add_option_group(group)

    (options, args) = parser.parse_args(sys.argv[1:])
    W.verbose = options.verbose
    W.tmppath = options.tmppath
    W.tmpbase = options.tmpbase
    W.keeptmp = options.keeptmp
    W.archive = options.archive
    W.archpatterns = options.archpatterns
    if options.stdout is not None:
        W.stdout = open(options.stdout,'w')
    if options.stderr is not None:
        W.stderr = open(options.stderr,'w')
    if options.stdin is not None:
        W.stdin = open(options.stdin,'r')
    W.setup()
    try:
        W.stagein()
        try:
            W.run()
        finally:
            W.stageout()
    finally:
        W.cleanup()

class Wrapper:
    def __init__(self):
        # User modifyable variables
        self.tmpbase = None
        self.tmppath = None
        self.keeptmp = False
        self.stdin = None
        self.stdout = None
        self.stderr = None
        self.cmdlines = []
        self.copyin = []
        self.copyout = []
        self.extract = []
        self.archive = None
        self.archpatterns = []
        self.verbose = False
        self.workdir = None
        self.env = os.environ.copy()

    def define_callback(self, option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
            self.env[ff[0]] = '='.join(ff[1:])
	else:
	    print >> sys.stderr, "Wrapper: Illegal argument to",opt," expected NAME=VALUE"
	    sys.exit(1)

    def infile_callback(self, option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
	    if len(ff) != 2:
                parser.error("Illegal argument to "+opt)
		print >> sys.stderr, "Wrapper: Illegal argument to",opt
		sys.exit(1)
	    self.infiles.append([ff[1],ff[0]])
	else:
	    self.infiles.append([value])

    def replace_callback(self, option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
	    ff = [ff[0],'='.join(ff[1:])]
	    if len(self.infiles) == 0:
		print >> sys.stderr, "Wrapper:",option,"given before any infile was specified"
		sys.exit(1)
	else:
	    print >> sys.stderr, "Wrapper: Illegal argument to",opt," expected NAME=VALUE"
	    sys.exit(1)
	self.infiles[-1].extend(ff)

    def getfile_callback(self, option, opt, value, parser):
	if '=' in value:
	    ff = value.split('=')
	    if len(ff) != 2:
		print >> sys.stderr, "Wrapper: Illegal argument to",opt
		sys.exit(1)
	    self.outfiles.append(ff)
	else:
	    self.outfiles.append([value])

    def command_callback(self, option, opt, value, parser):
	self.cmdlines.append(shlex.split(value))

    def modcmd_callback(self, option, opt, value, parser):
        "Accepts both --append and --prepend"
	if len(self.cmdlines) == 0:
		print >> sys.stderr, "Wrapper:",option,"given before any -c/--command was specified"
		sys.exit(1)
	if opt in ['-a','--append']:
	    self.cmdlines[-1].extend(shlex.split(value))
	else:
	    self.cmdlines[-1] = shlex.split(value) + self.cmdlines[-1]

    def setup(self):
        try:
            if self.tmppath is not None:
                if os.path.isdir(tmppath):
                    self.keeptmp = True
                else:
                    os.makedirs(self.tmppath)
            else:
                if self.tmpbase is None:
                    self.tmppath = tempfile.mkdtemp()
                else:
                    self.tmppath = tempfile.mkdtemp(dir=self.tmpbase)
            if self.verbose: print "Wrapper: Using temporary directory '%s'" % self.tmppath
            self.workdir = os.getcwd()
        except:
            print >> sys.stderr, "Wrapper: Cannot create temporary directory '%s', quitting." % self.tmppath
            raise

    def stagein(self):
        "Copy in files and extract archives. Stay in cwd"
        # Copy in files, stop if problem
        for f in self.copyin: 
            src = f[0]
            if len(f) % 2 == 0: # dst name given
                dst = os.path.join(self.tmppath,f[1])
                rep = f[2:]
            else:
                dst = os.path.join(self.tmppath,os.path.basename(src))
                rep = f[1:]
            rep = zip(rep[0::2],rep[1::2])
            try:
                self.copy_with_replacement(src,dst,rep)
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
        for a in self.extract:
            if self.verbose: print "Wrapper: Extracting archive",a
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
                t.extractall(path=self.tmppath)
            except:
                print >> sys.stderr, "Wrapper: Error extracting files from",a,"quitting."
                raise

    def copy_with_replacement(self,src,dst,rep=[]):
        try:
            if self.verbose: print "Wrapper: Copying '"+src+"' to '"+dst+"'"
            if len(rep) == 0:
                shutil.copy2(src,dst)
            else:
                if self.verbose:
                    for r in rep:
                        print "Wrapper: While replacing '"+r[0]+"' by '"+r[1]+"'"
                of = open(dst,'w')
                for l in open(src,'r').readlines():
                    for r in rep:
                        l = l.replace(r[0],r[1])
                    of.write(l)
                of.close()
        except IOError,e:
            print >> sys.stderr, "Wrapper: Error copying file '"+src+"'"
            raise    

    def stageout(self):
        # Create archive, easiest from tmppath
        os.chdir(self.tmppath)
        if self.archive is not None:
            try:
                self.archive = os.path.join(self.workdir,archive)
                if self.verbose: print "Wrapper: Creating archive",self.archive
                t = tarfile.open(self.archive,'w:gz')
                for p in self.archpatterns:
                    for f in glob.glob(p):
                        t.add(f)
                t.close()
            except:
               print >> sys.stderr, "Wrapper: Error creating tar archive"
               raise # Should try to copy out files as well
        # Copy back files
        os.chdir(self.workdir)
        def cpout(src,dst=None):
            if self.verbose: print "Wrapper: Copying back",src
            if os.path.isfile(src):
                if dst is None: dst='.'
                shutil.copy2(src,dst)
            elif os.path.isdir(src):
                if dst is None: dst=os.path.basename(src)
                if os.path.exists(dst): 
                    shutil.rmtree(dst)
                shutil.copytree(src,dst)
            elif self.verbose: # Just give a warning, don't panic
                print "Wrapper: When copying out, file",src,"not found."
        for f in self.copyout:
            try:
                if len(f) == 2:
                    cpout(os.path.join(self.tmppath,f[0]),f[1])
                else: # allow patterns in copyout if there is no = 
                    for ff in glob.glob(os.path.join(self.tmppath,f[0])):
                        cpout(ff)
            except:
                print >> sys.stderr,"Wrapper: Error copying out file",f
                raise # TODO: Must try to copy out all files

    def cleanup(self):
        if not self.keeptmp:
            shutil.rmtree(self.tmppath)
            if self.verbose: print "Wrapper: Temporary directory deleted"
        else:
            if self.verbose: print "Wrapper: Temporary directory kept:",tmppath

    def run(self):
        "Run cmdlines, stop if one fails."
        for cmd in self.cmdlines:
            pgid = None
            oldhand = []
            try:
                os.chdir(self.tmppath)
                cmd = "'"+"' '".join(cmd)+"'"
                def handler(signum,frame):
                    if pgid is not None: 
                        if signum == signal.SIGALRM: signum = signal.SIGTERM
                        os.killpg(pgid,signum)
                catch = [signal.SIGTERM, signal.SIGHUP, signal.SIGINT, signal.SIGALRM]
                for s in catch:
                    oldhand.append(signal.signal(s,handler))
                if self.verbose:
                    print "Wrapper: Executing",cmd
                p = subprocess.Popen(cmd,stdin=self.stdin,stdout=self.stdout,stderr=self.stderr,shell=True,
                                     preexec_fn=os.setsid,env=self.env)
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
                os.chdir(self.workdir)
        

wrap_dirac()
