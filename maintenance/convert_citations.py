#!/usr/bin/env python

"""
This script recursively replaces citations in reStructuredText files from this format:

    [Kullie2011]_

to the sphinx-contrib-bibtex compatible format:

    :cite:`Kullie2011`

Usage:

convert_citations.py /path/to/documentation/root /path/to/references/file.bib
"""

import os
import re
import sys
import logging
import argparse
import fileinput

from pybtex.database.input import bibtex

CITATION_REGEX = re.compile(r'\[([^\]]+)\]_')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source_path', help='specify the source path', type=str)
    parser.add_argument('bibtex_file', help='specify the bibtex file', type=str)
    args = parser.parse_args()

    source_path = os.path.abspath(args.source_path)
    if not os.path.isdir(source_path):
        parser.error('Invalid source path provided.')

    bibtex_file = os.path.abspath(args.bibtex_file)
    if not os.path.isfile(bibtex_file):
        parser.error('Invalid bibtex file provided.')

    logger = get_logger()
    citations = get_citations(bibtex_file)
    results = []

    recurse(source_path, results, citations, logger)

    logger.info('Found {0} results, replacing ...'.format(len(results)))

    for result in results:
        for line in fileinput.input(result['path'], inplace=True):
            print line.replace(result['replace_from'], result['replace_to']),


def get_logger():
    """
    Returns configured logger.
    """
    logger = logging.getLogger()
    formatter = logging.Formatter('%(message)s')

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.setLevel(logging.INFO)

    return logger


def get_citations(bibtex_file):
    """
    Returns citation keys of given BiBteX file.
    """
    parser = bibtex.Parser()
    data = parser.parse_file(bibtex_file)
    return data.entries.keys()


def recurse(path, results, citations, logger):
    """
    Recursively reads path and adds replacements to results list.
    """
    for item in os.listdir(path):
        item_path = os.path.join(path, item)

        if os.path.isdir(item_path):
            recurse(item_path, results, citations, logger)
        else:
            if item_path.endswith('.rst'):
                build_replacements(item_path, results, citations, logger)


def build_replacements(path, results, citations, logger):
    """
    Builds replacements for given path.
    """
    with open(path, 'r') as fp:
        contents = fp.read()

    for item in CITATION_REGEX.findall(contents):
        if item in citations:
            add_replace(path, results, item)
        else:
            logger.info('Citation "{0}" not found. Skipping ...'.format(item))


def add_replace(path, results, item):
    """
    Adds a new replacement for given path.
    """
    replace_from = '[{0}]_'.format(item)
    replace_to = ':cite:`{0}`'.format(item)

    results.append({
        'path': path,
        'item': item,
        'replace_from': replace_from,
        'replace_to': replace_to,
    })

if __name__ == '__main__':
    main()
