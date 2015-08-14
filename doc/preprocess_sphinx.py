#!/usr/bin/env python

"""
    Script to replace homecooked macros to get index and crossref labels
    'for free'. Bash version written by Radovan Bast, translated
    to Python by Ivan Hrasko. Further adaptation by Miro Ilias.

    Replaced macros in the rst-sources of the documentation are:
           keyword(SOMEKEYWORD)
           star(SOMEKEYWORD)
           starstar(SOMEKEYWORD)


"""

import re
import sys
import os
import fnmatch
import fileinput
import subprocess
import shutil
import string


def main(argv):

    KEYWORD_REGEX  = r'keyword\(([a-zA-Z !0-9_+-]+)\)'
    STARSTAR_REGEX = r'starstar\(([a-zA-Z 0-9_]+)\)'
    STAR_REGEX     = r'star\(([a-zA-Z 0-9_]+)\)'

    DIR_DOC   = os.path.abspath(argv[1])
    DIR_BUILD = os.path.abspath(argv[2])
    DIR_TEMP  = os.path.join(DIR_BUILD, 'html_temp')
    DIR_FINAL = os.path.join(DIR_BUILD, 'html')

    shutil.rmtree(DIR_TEMP, ignore_errors=True)
    shutil.copytree(DIR_DOC, DIR_TEMP)

    for root, dirs, files in os.walk(os.path.join(DIR_TEMP, 'manual')):
        for filename in fnmatch.filter(files, '*.rst'):

            # replace (in given **SECTION/*SECTION):
            #
            # keyword(RABOOF)
            # by:
            # .. index:: .RABOOF
            # .. SECTION_RABOOF:
            #
            # .RABOOF
            # -------
            section = ''
            for line in fileinput.input(os.path.join(root, filename), inplace=True):
                # the section is used to make the labels unique
                section_found = re.findall(STARSTAR_REGEX, line)
                if section_found == []:
                    section_found = re.findall(STAR_REGEX, line)
                if section_found != []:
                    section = '%s_' % string.replace(section_found[0], ' ', '_')
                print re.sub(KEYWORD_REGEX,
                            ('.. index:: .\g<1>\n.. _%s.\g<1>:\n\n.\g<1>\n'+str('-'*(line.find(')')-line.find('(')))) % section,line[0:len(line)-1])

            # replace:
            # starstar(RABOOF)
            # by:
            # .. index:: **RABOOF
            # .. _RABOOF:
            #
            # ==========
            # \*\*RABOOF
            # ==========
            for line in fileinput.input(os.path.join(root, filename), inplace=True):
                print re.sub(STARSTAR_REGEX,
                             ('.. index:: **\g<1>\n.. _**\g<1>:\n\n'+str('='*(line.find(')')-line.find('(')+3))+'\n\\*\\*\g<1>\n'+str('='*(line.find(')')-line.find('(')+3))),line[0:len(line)-1])

            # replace:
            # star(RABOOF)
            # by:
            # .. index:: **RABOOF
            # .. _RABOOF:
            #
            # ========
            # \*RABOOF
            # ========
            for line in fileinput.input(os.path.join(root, filename), inplace=True):
                print re.sub(STAR_REGEX,
                             ('.. index:: *\g<1>\n.. _*\g<1>:\n\n'+str('='*(line.find(')')-line.find('(')+1))+'\n\\*\g<1>\n'+str('='*(line.find(')')-line.find('(')+1))),line[0:len(line)-1])

    subprocess.call('sphinx-build -b html %s %s' % (DIR_TEMP, DIR_FINAL), shell=True)


if __name__ == '__main__':
    main(sys.argv)
