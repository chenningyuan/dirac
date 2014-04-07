# coding=utf-8

import os

from fabric.decorators import task
from fabric.contrib.files import exists
from fabric.contrib.console import confirm
from fabric.api import cd, env, prefix, run, local

env.repository = 'git@bitbucket.org:miroslav_ilias/my-dirac.git',
env.remote = 'origin'
env.branch = 'bubniak_dirac_devel'
env.user = 'jbubniak'
env.hosts = ['labs.fpv.umb.sk']
env.port = 22
env.root = '/home/jbubniak/my-dirac'
env.build_dir = 'build_bubniak'
env.build_path = '{0}/{1}'.format(env.root, env.build_dir)
env.html_path = '{0}/{1}/html'.format(env.root, env.build_dir)
env.public_html = '/home/jbubniak/public_html'


@task
def pull():
    with cd(env.root):
        branch = run('git rev-parse --abbrev-ref HEAD').strip()

        if branch != env.branch:
            run('git checkout -b {0} {1}/{2}'.format(env.branch, env.remote, env.branch))

        run('git pull {0} {1}'.format(env.remote, env.branch))


@task
def build():
    if exists(env.build_path):
        run('rm -rf {0}'.format(env.build_path))

    with cd(env.root):
        run('./setup {0}'.format(env.build_dir))


@task
def make_html():
    with cd(env.build_path):
        run('make html')


@task
def symlink_html():
    symlink_path = '{0}/html'.format(env.public_html)

    if not exists(symlink_path):
        run('ln -s {0} {1}'.format(env.html_path, env.public_html))


@task
def update():
    pull()
    build()
    make_html()
    symlink_html()
