# Universal wrapper, by Ulf Ekstrom 2011
import os, subprocess, shutil, glob, signal, optparse, sys, shlex

def _main_():
        usage = \
            """%prog [options] """
        parser =optparse.OptionParser(usage)

        group = optparse.OptionGroup(parser, 'General options')
        group.add_option('-v','--verbose', action='store_true', dest="verbose",
                         help='Print more information')
        group.add_option('-j','--jobname', type='string', default=None,dest="jobname",
                         help='Name of job, used as base for temporary directory and output files',metavar="JOBNAME")
        group.add_option('-t','--tmpbase', type='string', default=None,dest="tmpbase",
                         help='Base directory for temporary directory, default $TMPDIR',metavar="PATH")
        parser.add_option_group(group)

        group = optparse.OptionGroup(parser, 'Input options')
        infiles=[]
        def infilecb(option, opt, value, parser):
            if '=' in value:
                ff = value.split('=')
                if len(ff) != 2:
                    print "Illegal argument to",opt
                    sys.exit(1)
                infiles.append([ff[1],ff[0]])
            else:
                infiles.append([value])
        def replacecb(option, opt, value, parser):
            if '=' in value:
                ff = value.split('=')
                ff = [ff[0],'='.join(ff[1:])]
                if len(infiles) == 0:
                    print option,"given before any -i/--infile was specified"
                    sys.exit(1)
            else:
                print "Illegal argument to",opt," expected NAME=VALUE"
                sys.exit(1)
            infiles[-1].extend(ff)
        group.add_option('-i','--infile', type='string', action='callback',callback=infilecb,
                         help='Copy file to TMPDIR/NAME',metavar='[NAME=]FILE')
        group.add_option('-r','--replace', type='string', action='callback',callback=replacecb,
                         help='Replace NAME by VALUE in the previous --infile.',metavar="NAME=VALUE")
        parser.add_option_group(group)

        group = optparse.OptionGroup(parser, 'Command options')
        commands=[]
        def commandcb(option, opt, value, parser):
            commands.append(shlex.split(value))
        def pendcb(option, opt, value, parser):
            if len(commands) == 0:
                    print option,"given before any -c/--command was specified"
                    sys.exit(1)
            if opt in ['-a','--append']:
                commands[-1].extend(shlex.split(value))
            else:
                commands[-1] = shlex.split(value) + commands[-1]
        group.add_option('-c','--command', type='string', action='callback',callback=commandcb,
                         help='Add a new command to execute',metavar='COMMAND')
        group.add_option('-a','--append', type='string', action='callback',callback=pendcb,
                         help='add ARGS last on the previous command line given by -c/--command',metavar='ARGS')
        group.add_option('-p','--prepend', type='string', action='callback',callback=pendcb,
                         help='add ARGS first on the previous command line given by -c/--command',metavar='ARGS')
        parser.add_option_group(group)

        (options, args) = parser.parse_args(sys.argv[1:])
        
	if options.verbose:
            print 'DEBUG: infiles',infiles
	    print 'DEBUG: commands',commands
        if options.jobname is None:
            print "No job name (-j) given, quitting"
            sys.exit(1)
        if options.tmpbase is None:
            for d in ['TMPDIR','TMP']:
                if os.environ.has_key(d):
                    options.tmpbase = os.environ[d]
                    break
        tmppath = os.path.join(options.tmpbase,options.jobname+'.'+str(os.getpid()))
        try:
            run_wrapped(tmppath,commands,copyin=infiles,verbose=options.verbose)
        except:
            print "Something went wrong, re-raising Python error.."
            raise

def run_wrapped(tmppath, cmdlines, stdin=None, stdout=None, stderr=None,
                copyin=[], copyout=[], verbose=False, 
                delete_tmp=True):
    """Create tmpdir, cd there and execute cmdlines, a list of arguments
    to Popen.
    cmdlines is a list of command lines in list form, ex: cmdlines=[['ls','-l'],[..]]
    copyin is a list of files to copy, in the form
    copyin=[(frompath,localname,rep1,with1,rep2,with2,..)] where rep1 will be 
    replaced by with1 during the copy.
    localname can be omitted, then the basename of frompath will be used.
    copyout is similar but no replacement is done
    """
    # Create temporary directory
    workdir = os.getcwd()
    if os.path.isdir(tmppath):
        old_tmpdir = True
        if verbose: print "Temporary directory '%s' already there" % tmppath
    else:
        try:
            os.makedirs(tmppath)
        except:
            print "Cannot create temporary directory '%s', quitting." % tmppath
            raise
        old_tmpdir = False        
    try: # Need to delete the tmp dir after this        
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
                print "File copy failed, quitting."
                raise
        try: # Run commands, stop if one fails
            os.chdir(tmppath)
            for cmd in cmdlines:
                pgid = None
                oldhand = []
                try:
                    cmd = "'"+"' '".join(cmd)+"'"
                    def handler(signum,frame):
                        if pgid is not None: os.killpg(pgid,signum)
                    catch = [signal.SIGTERM, signal.SIGHUP, signal.SIGINT, signal.SIGALRM]
                    for s in catch:
                        oldhand.append(signal.signal(s,handler))
                    if verbose:
                        print "Executing",cmd
                    p = subprocess.Popen(cmd,stdin=stdin,stdout=stdout,stderr=stderr,shell=True,
                                         preexec_fn=os.setsid)
                    pgid = p.pid
                    p.wait() # OSError will be raised here if killed
                    if p.returncode != 0:
                        raise subprocess.CalledProcessError(p.returncode,cmd)
                except subprocess.CalledProcessError, e:
                    print "Command "+cmd+" failed with returncode",e.returncode
                    raise e
                except OSError, e: # Typically signal
                    print "Got OSError %i while executing" % e.errno,cmd
                    raise e
                finally: # restore signal handlers
                    if len(oldhand) > 0:
                        for s,h in zip(catch,oldhand):
                            signal.signal(s,h)
        finally: # Copy back files, always do this even if a command fails.
            os.chdir(workdir)
            for f in copyout:
                src = os.path.join(tmppath,f[0])
                if os.path.isfile(src): # It's ok if the file is not there
                    if len(f) == 2: # Destination path given
                        dst = f[1]
                    else:
                        dst = os.path.basename(src)
                    try:
                        copy_with_replacement(src,dst,[],verbose=verbose)
                    except:
                        print "File copy out of '%s' to '%s' failed, continuing" % (src,dst)
                elif verbose: # Just give a warning
                    print "When copying out, file",src,"not found."

            # Archive TODO

    finally:  # Delete temporary directory
        if delete_tmp and not old_tmpdir: # If tmpdir was there already don't delete it
            shutil.rmtree(tmppath)
            if verbose: print "Temporary directory deleted"
        else:
            if verbose: print "Temporary directory kept:",tmppath

def copy_with_replacement(src,dst,rep=[],verbose=False):
    try:
        if verbose: print "Copying '"+src+"' to '"+dst+"'"
        if len(rep) == 0:
            shutil.copyfile(src,dst)
        else:
            if verbose:
                for r in rep:
                    print " Replacing '"+r[0]+"' by '"+r[1]+"'"
            of = open(dst,'w')
            for l in open(src,'r').readlines():
                for r in rep:
                    l = l.replace(r[0],r[1])
                of.write(l)
            of.close()
    except IOError,e:
        print "Error copying file '"+src+"'"
        raise        

_main_()

