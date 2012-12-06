#!/usr/bin/env python
import sys, urllib2, subprocess, os

# Update script 2012 by Ulf Ekstrom

# The url for the next update is read from this file
url_location = 'utils/update-url'

def get_patch(url):
    # print "Attempting to download patch from",url
    try:
        u = urllib2.urlopen(url,timeout=10)
    except urllib2.HTTPError:
        print "No updates found on server."
        return None
    except urllib2.URLError:
        sys.stderr.write("Error downloading patch from %s\n" % url)
        return None
    print "Downloading patch:",
    try:
        patch = u.read()
    except:
        sys.stderr.write("Error downloading patch from %s\n" % url)
        return None
    print "OK"
    return patch

def has_patchprogram():
    with open(os.devnull, "w") as fnull:
        res = subprocess.call(['patch','-v'],stdout = fnull, stderr = fnull)
    if res == 0:
        return True
    else:
        return False

def test_or_apply_patch(patch,test=True):
    with open(os.devnull, "w") as fnull:
        args = ['patch','-p1']
        if test:
            args.append('--dry-run')
        p = subprocess.Popen(args,stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.communicate(patch)
        res = p.returncode
    if res == 0:
        return True
    else:
        f = open('update-error.log','w')
        f.write(out[0])
        f.write(out[1])
        return False

def main():
    try:
        f = open(url_location,'r')
        url = f.read().strip()
    except:
        sys.stderr.write("Cannot find the file `%s', running update from wrong directory?\n" % url_location)
        sys.exit(-1)
    p = get_patch(url)
    if p:
        if not has_patchprogram():
            sys.stderr.write("The `patch' program was not found, cannot apply updates.\n")
            sys.exit(-1)
        if test_or_apply_patch(p,test=True):
            print "Verified patch, now applying",
            res = test_or_apply_patch(p,test=False)
            if res:
                print "Update applied successfully, please rerun setup.py and recompile Dirac"
                try:
                    f = open(url_location,'r')
                    newurl = f.read().strip()
                except:
                    sys.stderr.write("After update: Cannot find the file `%s', please report this to the developers.\n" % url_location)
                    sys.exit(-1)
                if newurl == url:
                    sys.stderr.write("After update: Idiotic patch, new update url same as the old. Please report this to the developers.\n")
                    sys.exit(-1)
            else:
                sys.stderr.write("Something went wrong applying the patch even after verification, see update-error.log.\n")
                sys.stderr.write("Please extract the official distribution archive and update directly\nfrom that.\n")
                sys.stderr.write("The source tree may be left in an inconsistent and non-compiling state.\n")
                sys.exit(-1)
        else:
            sys.stderr.write("Patch cannot be applied, please find error log in update-error.log.\n")
            sys.stderr.write("This may be due to your modifications of Dirac or a developer error.\n")
            sys.stderr.write("Please extract the official distribution archive and update directly\nfrom that.\n")
            sys.stderr.write("The error was detected before any files where modified.\n")
            sys.exit(-1)
    else:
        sys.exit(0)

print 'Dirac update script, attempting to locate and download updates..'
print
while True:
    main()
    print "Checking for more updates.."
