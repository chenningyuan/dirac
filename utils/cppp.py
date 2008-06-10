#! /usr/bin/env python

import sys, re, os

class FriendlyFile:
    def __init__(self, path, mode):
        self.f = open(path,mode)
        self.line_nr = 0
        self.path = path
    def readline(self):
        l = self.f.readline()
        if len(l) > 0:
            self.line_nr += 1
        # sys.stderr.write('Read line %i: %s\n' % (self.line_nr,l))
        return l
    def locstring(self):
        return "%s:%i" % (self.path,self.line_nr)
    

re_ifdef = re.compile('^\s*#\s*ifdef\s+(\w*)')
re_ifndef = re.compile('^\s*#\s*ifndef\s+(\w*)')
re_else = re.compile('^\s*#\s*else')
re_endif = re.compile('^#\s*endif')
re_if = re.compile('^\s*#\s*if\s+(\S*)')
re_elif = re.compile('^\s*#\s*elif\s+(\S*)')

defined = []
undefined = []


def process_chunk(infile, dst, do_echo = 1, cond_active = 0):
    line = infile.readline()
    while len(line) > 0:
        # Warn in case of unsupported lines
        # TODO: Must match whole token instead of string.find()
#         if re_if.match(line) or re_elif.match(line):
#             found = 0
#             for n in defined:
#                 if line.find(n) > -1:
#                     found = 1
#                     break
#             for n in undefined:
#                 if line.find(n) > -1:
#                     found = 1
#                     break
#             if found:
#                 sys.stderr.write("Warning at %s\n" % infile.locstring())
#                 sys.stderr.write("The line '%s' seem to contain one of your active symbols. This is not supported.\n" % line.strip())  
        # Actual work

        if re_ifdef.match(line):
            m = re_ifdef.match(line)
            name = m.group(1)
            if name in defined:
                process_chunk(infile,dst,do_echo, cond_active = 1)
            elif name in undefined:
                process_chunk(infile,dst,do_echo = 0, cond_active = 1)
            else:
                if do_echo: dst.write(line)
                process_chunk(infile,dst,do_echo, cond_active = 0)
        elif re_ifndef.match(line):
            m = re_ifndef.match(line)
            name = m.group(1)
            if name in undefined:
                process_chunk(infile,dst,do_echo, cond_active = 1)
            elif name in defined:
                process_chunk(infile,dst,do_echo = 0, cond_active = 1)
            else:
                if do_echo: dst.write(line)
                process_chunk(infile,dst,do_echo, cond_active = 0)
        elif re_if.match(line):
            if do_echo: dst.write(line)
            process_chunk(infile,dst,do_echo, cond_active = 0)
        elif re_else.match(line):
            if cond_active:
                do_echo = 1 - do_echo
            elif do_echo:
                dst.write(line)
        elif re_endif.match(line):
            if (not cond_active) and do_echo:
                dst.write(line)
            return
        elif do_echo:
            dst.write(line)
                      
        line = infile.readline()

usage = """Usage: cppp [-DMACRO ..] [-UMACRO ...] FILE [FILE ..]
For each FILE perform a partial C preprocessing of #ifdef and #ifndef
directives, completely removing sections that cannot be reached with
the given macros defined (-D) or undefined (-U).

For example, the file

#ifdef A
So A was defined
#ifdef C
#ifndef B
Only if B is not defined
#endif
#endif
#else
A was not defined
#endif

With the arguments -DA -UB will be transformed to

So A was defined
#ifdef C
Only if B is not defined
#endif


NOTE: The program does not process #if statements, although they are reproduced
correctly in the preprocessed files. You may have to manually convert such statements
to chained #ifdef's.

NOTE2: The files are modified IN PLACE! Run this on a copy of the source tree.

By U. Ekstrom, <uekstrom@gmail.com> 2008

"""


def main():
    filenames = []
    target = filenames
    if len(sys.argv) == 1:
        print usage
        return    
    for arg in sys.argv[1:]:
        arg = arg.strip()
        if arg.find('-D') == 0:
            arg = arg[2:]
            target = defined
        elif arg.find('-U') == 0:
            arg = arg[2:]
            target = undefined
        if len(arg) > 0:
            target.append(arg)
            if target in [defined, undefined]:
                if '=' in arg:
                    print 'Error: Giving a macro a value (%s) is not supported. Quitting!' % arg
                    return -1
        target = filenames
    for fn in filenames:
        f = FriendlyFile(fn,'r')
        fout = open(fn+'.cppp~','w')
        process_chunk(f,fout)
        fout.close()
        os.rename(fn+'.cppp~',fn)

main()
