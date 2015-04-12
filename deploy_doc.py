#!/usr/bin/env python

import os
import json
import shutil
import argparse
import datetime
import subprocess


class Syncer(object):
    def __init__(self, root, user, host, port=22):
        self.root = root
        self.user = user
        self.host = host
        self.port = port

    def sync(self, src, dest):
        """
        Syncs the local directory with the remote one.
        """
        dest = os.path.join(self.root, dest)
        args = {
            'user': self.user,
            'host': self.host,
            'port': self.port,
            'src': src,
            'dest': dest,
        }
        cmd = 'rsync -arzvP -e "ssh -o StrictHostKeyChecking=no -p {port}" ' \
              '--rsync-path="mkdir -p {dest} && rsync" {src} {user}@{host}:{dest}'.format(**args)
        return subprocess.call(cmd, shell=True)


class Deployment(object):
    name = None
    log_name = 'build_log.txt'
    metadata_name = 'build_metadata.json'

    def __init__(self, build_dir, syncer):
        self.build_dir = build_dir
        self.syncer = syncer

    def deploy(self):
        """
        Deploys all of the builds.
        """
        branch = self.get_branch()
        revision = self.get_revision()

        self.deploy_sphinx(branch=branch, revision=revision)
        self.deploy_slides(branch=branch, revision=revision)
        self.deploy_doxygen(branch=branch, revision=revision)

    def deploy_sphinx(self, branch, revision):
        """
        Deploys the Sphinx build.
        """
        name = 'sphinx'
        src = self.get_module_src_path(name)
        dest = self.get_deploy_dest(module=name, branch=branch, revision=revision)
        log = self.write_log(module=name)

        self.write_metadata(src, module=name, branch=branch, revision=revision, log=log)

        return self.syncer.sync(src=src, dest=dest)

    def deploy_slides(self, branch, revision):
        """
        Deploys the slides build.
        """
        name = 'slides'
        src = self.get_module_src_path(name)
        dest = self.get_deploy_dest(module=name, branch=branch, revision=revision)
        log = self.write_log(module=name)

        self.write_metadata(src, module=name, branch=branch, revision=revision, log=log)

        return self.syncer.sync(src=src, dest=dest)

    def deploy_doxygen(self, branch, revision):
        """
        Deploys the doxygen build.
        """
        name = 'doxygen'
        src = self.get_module_src_path(name)
        dest = self.get_deploy_dest(module=name, branch=branch, revision=revision)
        log = self.write_log(module=name)

        self.write_metadata(src, module=name, branch=branch, revision=revision, log=log)

        return self.syncer.sync(src=src, dest=dest)

    def get_module_src_path(self, module):
        """
        Returns path to the builded HTML root directory for given module.
        """
        mapping = {
            'sphinx': ['html'],
            'slides': ['html_slides'],
            'doxygen': ['_doxygen', 'html'],
        }
        return os.path.join(self.get_build_path(), *mapping[module])

    def get_log_name(self, module):
        """
        Returns the expected build log filename for given module.
        """
        mapping = {
            'sphinx': 'sphinx-log.txt',
            'slides': 'slides-log.txt',
            'doxygen': 'doxygen-log.txt',
        }
        return mapping[module]

    def write_log(self, module):
        """
        Writes the build log of given module to the module build directory.
        """
        log_name = self.get_log_name(module)
        log_path = os.path.join(self.get_build_path(), log_name)

        if os.path.exists(log_path):
            dest = os.path.join(self.get_module_src_path(module), self.log_name)
            shutil.copy(src=log_path, dst=dest)
            return self.log_name

        return None

    def get_deploy_dest(self, module, branch, revision):
        """
        Returns path to the remote deploy destination. It will be joined with the
        syncer root path later.

        Example: semaphoreci/sphinx/master/6das54d65as4d65as465d4as6d4a5s
        """
        return os.path.join(self.name, module, branch, revision)

    def write_metadata(self, src, module, branch, revision, log):
        """
        Writes the metadata of the build in the JSON format.
        """
        metadata = self.build_metadata(module, branch, revision, log)
        metadata_path = os.path.join(src, self.metadata_name)
        with open(metadata_path, 'w') as fp:
            json.dump(obj=metadata, fp=fp, indent=4)
        return metadata

    def build_metadata(self, module, branch, revision, log):
        """
        Returns build metadata as a python dictionary.
        """
        return {
            'ci': self.name,
            'module': module,
            'branch': branch,
            'revision': revision,
            'message': self.get_commit_message(),
            'date': datetime.datetime.now().isoformat(),
            'log': log,
        }

    def get_build_path(self):
        """
        Returns path to the DIRAC build directory.
        """
        return os.path.join(self.get_root_path(), self.build_dir)

    def get_root_path(self):
        """
        Returns path to the root directory - git repository root.
        """
        raise NotImplementedError

    def get_revision(self):
        """
        Returns the current git revision.
        """
        raise NotImplementedError

    def get_branch(self):
        """
        Returns the current git branch.
        """
        raise NotImplementedError

    def get_commit_message(self):
        """
        Returns the git comment of the current revision.
        """
        args = {
            'root': self.get_root_path(),
            'revision': self.get_revision(),
        }
        cmd = 'cd {root} && git log --format=%B -n 1 {revision}'.format(**args)
        data = subprocess.check_output(cmd, shell=True)
        return data.strip() if data else ''


class Semaphore(Deployment):
    name = 'semaphoreci'

    def get_root_path(self):
        return os.environ['SEMAPHORE_PROJECT_DIR']

    def get_branch(self):
        return os.environ['BRANCH_NAME']

    def get_revision(self):
        return os.environ['REVISION']


def main():
    ci_mapping = {
        'semaphoreci': Semaphore,
    }

    parser = argparse.ArgumentParser()
    parser.add_argument('ci', help='specify the CI name', type=str)
    parser.add_argument('build_dir', help='specify the build dir name', type=str)
    parser.add_argument('--root', help='specify the remote root path', type=str)
    parser.add_argument('--user', help='specify the remote user name', type=str)
    parser.add_argument('--host', help='specify the remote host name', type=str)
    parser.add_argument('--port', help='specify the remote port', type=int)
    args = parser.parse_args()

    if not args.ci or args.ci not in ci_mapping:
        parser.error('Unsupported CI service {0}'.format(args.ci))

    if not args.build_dir:
        parser.error('Build dir is not valid.')

    deploy(ci_mapping, args)


def deploy(mapping, args):
    syncer = Syncer(root=args.root, user=args.user, host=args.host, port=args.port)

    deployer = mapping[args.ci](build_dir=args.build_dir, syncer=syncer)
    deployer.deploy()


if __name__ == '__main__':
    main()
