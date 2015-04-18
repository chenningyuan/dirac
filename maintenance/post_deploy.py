import os
import json
import datetime
import shutil
import itertools
import logging
import logging.handlers

from jinja2 import Template

# keep only limited number of builds per CI service
CI_BUILD_LIMIT = 5

# where are deployments stored
DEPLOY_PATH = '/var/www/html/sphinx/doc-deployments'

# path to generated index file
DATA_FILE = '/var/www/html/sphinx/doc-deployments/index.html'

# repo urls
BRANCH_URL = 'https://bitbucket.org/miroslav_ilias/my-dirac/commits/branch/{0}'
REVISION_URL = 'https://bitbucket.org/miroslav_ilias/my-dirac/commits/{0}'

# the main index template
TEMPLATE = """
    <!DOCTYPE html>
    <html>
        <head>
            <meta charset="utf-8">
            <meta http-equiv="X-UA-Compatible" content="IE=edge">
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <title>DIRAC Deployments</title>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css">
        </head>
        <body>
            <div class="container-fluid">
                <div class="row">
                    <div class="col-md-12">
                        <div class="page-header">
                            <h1>Deployments <small>DIRAC documentation</small></h1>
                        </div>
                        {% for ci, rows in data.items() %}
                            <div class="panel panel-default">
                                <div class="panel-heading">
                                    <h3 class="panel-title">{{ ci }}</h3>
                                </div>
                                <div class="panel-body">
                                    <div class="table-responsive">
                                        <table class="table table-striped table-condensed table-hover">
                                            <thead>
                                                <tr>
                                                    <th>Date</th>
                                                    <th>Branch</th>
                                                    <th>Revision</th>
                                                    <th>Message</th>
                                                    <th colspan="3">Modules</th>
                                                </tr>
                                            </thead>
                                            <tbody>
                                                {% for row in rows[:build_limit] %}
                                                    <tr>
                                                        <td>{{ row['rev']['metadata']['date'] }}</td>
                                                        <td><a href="{{ row['branch']['link'] }}">{{ row['branch']['name'] }}</a></td>
                                                        <td><a href="{{ row['rev']['link'] }}">{{ row['rev']['name']|truncate(6, True) }}</a></td>
                                                        <td>{{ row['rev']['metadata']['message']|truncate(60) }}</td>
                                                        {% for group in row['groups'] %}
                                                            <td>
                                                                <strong>{{ group['module']|capitalize }}</strong>
                                                                <span>
                                                                    {% if group['html_exists'] %}
                                                                        <a href="/deploy/{{ ci }}/{{ group['module'] }}/{{ row['branch']['name'] }}/{{ row['rev']['name'] }}/" class="btn btn-default btn-xs">HTML</a>
                                                                    {% endif %}
                                                                </span>
                                                                <span>
                                                                    <a href="/deploy/{{ ci }}/{{ group['module'] }}/{{ row['branch']['name'] }}/{{ row['rev']['name'] }}/{{ group['log'] }}" class="btn btn-default btn-xs">LOG</a>
                                                                </span>
                                                            </td>
                                                        {% endfor %}
                                                    </tr>
                                                {% endfor %}
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            </div>
                        {% endfor %}
                    </div>
                </div>
            </div>
        </body>
    </html>
"""


def main():
    logger = get_logger('deployer.txt')

    if not os.path.exists(DEPLOY_PATH):
        logger.info('Deploy path not found. Exiting ...')
        return

    logger.info('Starting ...')

    data = get_data()
    data = flatten_data(data)
    data = order_data(data)
    data = group_data(data)

    logger.info('Writing new data ...')

    write_data(data)

    logger.info('Cleaning old builds ...')

    clean_data(data, logger)

    logger.info('Finishing ...')


def get_logger(filename):
    """
    Returns configured logger to log to given filename.
    """
    logger = logging.getLogger()
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    # limit to 1MB, keep 5 last records
    file_handler = logging.handlers.RotatingFileHandler(filename, maxBytes=1048576, backupCount=5)
    file_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)

    return logger


def write_data(data):
    """
    Renders data using HTML template and writes it to file.
    """
    with open(DATA_FILE, 'w') as fp:
        template = get_template(data)
        fp.write(template)


def clean_data(data, logger):
    """
    Removes the old builds.
    """
    removed = 0

    for ci, rows in data.items():
        build_count = len(rows)

        if build_count > CI_BUILD_LIMIT:
            builds_to_remove = rows[CI_BUILD_LIMIT:]
            removed += remove_builds(ci, builds_to_remove, logger)

    if removed:
        logging.info('Removed {0} builds ...'.format(removed))
    else:
        logging.info('No need to remove any builds ...')


def remove_builds(ci, builds, logger):
    """
    Removes given builds from filesystem.
    """
    removed = 0

    for build in builds:
        for group in build['groups']:
            metadata_path = get_metadata_path(group['path'])

            # double check if metadata exists in build path
            if os.path.exists(metadata_path):
                logging.info('Removing build "{0}" ...'.format(group['path']))

                # double check if exists
                if os.path.exists(group['path']):
                    shutil.rmtree(group['path'])

    return removed


def get_template(data):
    """
    Returns rendered Jinja2 template with given data as context.
    """
    template = Template(TEMPLATE.strip())
    return template.render(data=data, build_limit=CI_BUILD_LIMIT)


def flatten_data(data):
    """
    Returns flattened data list.
    """
    flatten = []
    for ci, ci_data in data.items():
        for mod, mod_data in ci_data['modules'].items():
            for br, br_data in mod_data['branches'].items():
                for rev, rev_data in br_data['revisions'].items():
                    row = build_flat_row(ci, ci_data, mod, mod_data, br, br_data, rev, rev_data)
                    flatten.append(row)
    return flatten


def build_flat_row(ci, ci_data, mod, mod_data, br, br_data, rev, rev_data):
    """
    Builds flattend row from various data.
    """
    return {
        'ci': {
            'name': ci,
            'path': ci_data['path'],
        },
        'module': {
            'name': mod,
            'path': mod_data['path'],
        },
        'branch':  {
            'name': br,
            'path': br_data['path'],
            'link': br_data['link'],
        },
        'rev':  {
            'name': rev,
            'path': rev_data['path'],
            'link': rev_data['link'],
            'metadata': rev_data['metadata'],
        },
    }


def order_data(data):
    """
    Returns data ordered by build date.
    """
    return sorted(data, key=lambda row: row['rev']['metadata']['date'], reverse=True)


def group_data(data):
    """
    Returns data grouped by CI name.
    """
    d = {}
    for key, group in itertools.groupby(data, key=lambda row: row['ci']['name']):
        if key not in d:
            d[key] = []
        for r in group_rev_data([g for g in group]):
            d[key].append(r)
    return d


def group_rev_data(data):
    """
    Groups rows by revision.
    """
    l = []
    for key, group in itertools.groupby(data, key=lambda row: row['rev']['name']):
        groups = list(group)
        row = groups[0]
        row['groups'] = [build_row_group(g) for g in groups]
        l.append(row)
    return l


def build_row_group(row):
    """
    Returns row group.
    """
    metadata = row['rev']['metadata']
    metadata['date'] = metadata['date'].strftime('%Y-%m-%d %H:%M:%S')
    return {
        'module': row['module']['name'],
        'path': row['rev']['path'],
        'log': metadata['log'],
        'html_exists': metadata['html_exists'],
    }


def get_data():
    """
    Reads all of the build data from filesystem.
    """
    data = {}
    cis = os.listdir(DEPLOY_PATH)

    for ci in cis:
        ci_path = os.path.join(DEPLOY_PATH, ci)

        if not os.path.isdir(ci_path):
            continue

        data[ci] = {
            'path': ci_path,
            'modules': {},
        }

        for module in os.listdir(ci_path):
            module_path = os.path.join(ci_path, module)

            if not os.path.isdir(module_path):
                continue

            data[ci]['modules'][module] = get_module_data(module, module_path)

            for branch in os.listdir(module_path):
                branch_path = os.path.join(module_path, branch)

                if not os.path.isdir(branch_path):
                    continue

                data[ci]['modules'][module]['branches'][branch] = \
                    get_branch_data(branch, branch_path)

                for revision in os.listdir(branch_path):
                    revision_path = os.path.join(branch_path, revision)

                    if not os.path.isdir(revision_path):
                        continue

                    metadata = get_metadata(module, revision_path)

                    # skip build with no metadata
                    if not metadata:
                        continue

                    data[ci]['modules'][module]['branches'][branch]['revisions'][revision] = \
                        get_revision_data(revision, revision_path, metadata)

    return data


def get_module_data(module, module_path):
    """
    Returns module data.
    """
    return {
        'path': module_path,
        'branches': {},
    }


def get_branch_data(branch, branch_path):
    """
    Returns branch data.
    """
    return {
        'path': branch_path,
        'link': get_branch_link(branch),
        'revisions': {},
    }



def get_revision_data(revision, revision_path, metadata):
    """
    Returns revision data.
    """
    return {
        'path': revision_path,
        'link': get_revision_link(revision),
        'metadata': metadata,
    }


def get_metadata(module, revision_path):
    """
    Reads the build metadata from JSON file.
    """
    metadata_path = get_metadata_path(revision_path)
    if os.path.exists(metadata_path):
        with open(metadata_path, 'r') as fp:
            data = json.load(fp)
            data['date'] = convert_date(data['date'])
            data['html_exists'] = html_exists(revision_path)
            return data
    return None


def get_metadata_path(revision_path):
    """
    Returns path to the metadata file of given revision path.
    """
    return os.path.join(revision_path, 'build_metadata.json')


def convert_date(date_str):
    """
    Converts date string to date using iso format.
    """
    return datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%f')


def html_exists(revision_path):
    """
    Checks whether the HTML build was successfull and produced index.html file.
    """
    return os.path.exists(os.path.join(revision_path, 'index.html'))


def get_branch_link(branch):
    """
    Returns link to branch in VCS project website.
    """
    return BRANCH_URL.format(branch)


def get_revision_link(revision):
    """
    Returns link to revision in VCS project website.
    """
    return REVISION_URL.format(revision)

if __name__ == '__main__':
    main()
