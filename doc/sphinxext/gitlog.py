# coding=utf-8

import os
import re
import subprocess

from datetime import datetime

DEFAULT_DATE_FORMAT = '%d.%m.%Y %H:%M:%S'
LOG_FORMAT = 'name=%an email=%ae date=%at sha=%H'
LOG_REGEX = re.compile(u'^name=(?P<name>[\w\s]+) email=(?P<email>[\w@\.\s\-]+) date=(?P<date>\d{10}) sha=(?P<sha>\w{40})$', re.UNICODE)


def git_installed():
    process = subprocess.Popen('git --version', shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.communicate()
    return process.returncode == 0


def get_git_command(page_path, log_format):
    return 'git log -1 --format="{0}" -- {1}'.format(log_format, page_path)


def get_page_path(cfg, pagename):
    return os.path.join(os.path.abspath(cfg.gitlog_source_path), '{0}.rst'.format(pagename))


def get_git_log(cmd):
    return subprocess.check_output(cmd, shell=True)


def parse_log_info(info):
    match = LOG_REGEX.search(info.strip())
    if not match:
        print info
    return match.groupdict() if match else None


def get_git_info(cfg, pagename):
    page_path = get_page_path(cfg, pagename)
    cmd = get_git_command(page_path, LOG_FORMAT)
    info = get_git_log(cmd)
    return parse_log_info(info) if info else None


def format_date_string(date_str, date_format):
    try:
        dt = datetime.fromtimestamp(float(date_str))
    except ValueError:
        return None
    return dt.strftime(date_format)


def fill_context(cfg, context, info):
    context['git_author_name'] = info.get('name', '')
    context['git_author_email'] = info.get('email', '')
    context['git_commit_sha'] = info.get('sha', '')
    context['git_commit_date'] = format_date_string(info.get('date', ''), cfg.gitlog_date_format)


def page_context_handler(app, pagename, templatename, context, doctree):
    if git_installed():
        info = get_git_info(app.config, pagename)

        if info:
            fill_context(app.config, context, info)


def setup(app):
    app.add_config_value('gitlog_source_path', DEFAULT_DATE_FORMAT, 'html')
    app.add_config_value('gitlog_date_format', DEFAULT_DATE_FORMAT, 'html')
    app.connect('html-page-context', page_context_handler)