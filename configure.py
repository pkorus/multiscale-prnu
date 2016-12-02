#!/usr/bin/python
import commands
import sys
import os
import urllib2
import math

deps = {
    'ugm': 'http://www.cs.ubc.ca/~schmidtm/Software/UGM_2011.zip',
    'maxflow': 'http://vision.csd.uwo.ca/code/maxflow-v3.01.zip',
    'dde-prnu': 'http://dde.binghamton.edu/download/camera_fingerprint/CameraFingerprint_1.0.zip'
}

post = {
    'ugm': 'mv 3rd-party/ugm/UGM/* 3rd-party/ugm && rm -r 3rd-party/ugm/UGM',
    'maxflow': 'mv 3rd-party/maxflow 3rd-party/max-flow/maxflow-v3.0',
    'dde-prnu': 'mv 3rd-party/dde-prnu/CameraFingerprint/* 3rd-party/dde-prnu && rm -r 3rd-party/dde-prnu/CameraFingerprint',
    'tifs-2016-maps' : 'mv data/tifs-2016-maps/data-tifs-2016-maps/* data/tifs-2016-maps && rm -r data/tifs-2016-maps/data-tifs-2016-maps',
    'images' : 'mv data/images/data-images/* data/images && rm -r data/images/data-images && mv data/images/camera_models/* data/camera_models && rm -r data/images/camera_models'
}

sizes = {'UGM_2011.zip': 477992, 'CameraFingerprint_1.0.zip': 2113683, 'maxflow-v3.01.zip': 15006,
         'data-tifs-2016-maps.zip': 147930670, 'realistic-tampering-dataset.zip': 1697070982}


def match_command(input):
    patterns = ['clean', 'dependencies', 'data:maps', 'data:images']
    matching = [p for p in patterns if p.startswith(input)]
    if len(matching) == 1: return matching[0]
    return ''


def print_usage():
    print('usage: ./configure.py (clean|dependencies|data:maps|data:images)\n')
    print('  clean        - cleans up 3-rd party dependencies')
    print('  dependencies - download and install 3rd-party dependencies')
    print('  data:maps    - downloads tampering localization maps (140 MB)')
    print('  data:images  - downloads full size images (1.3 GB)')
    print('\nPartial commands are supported, e.g., ./configure dep -> ./configure dependencies')


def sh(cmd):
    # print '#> %s' % cmd
    return commands.getstatusoutput(cmd)


def unzip(filename, dirname):
    return sh('unzip %s -d %s' % (filename, dirname))


def download(url, filename, dirname='3rd-party'):
    u = urllib2.urlopen(url)
    f = open('%s/%s' % (dirname, filename), 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print('Downloading: %s (%sb) [' % (filename, file_size)),

    file_size_dl = 0
    block_sz = 8192
    progress = 0

    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)
        if math.floor(file_size_dl * 100. / file_size / 5) > progress:
            sys.stdout.write('#')
            sys.stdout.flush()
            progress = math.floor(file_size_dl * 100. / file_size / 5)

    print(']')
    f.close()


def isfileok(filename, dirname):

    if not os.path.exists('%s/%s' % (dirname,filename)):
        return False

    file_size = os.path.getsize('%s/%s' % (dirname,filename))
    if file_size != sizes[filename]:
        print('Cached copy of %s is invalid, cleaning up...' % filename)
        os.remove('%s/%s' % (dirname,filename))
        return False

    return True


def fetch_dependency(dep, url, cache_dir):
    dir_name = '%s/%s' % (cache_dir, dep)
    filename = os.path.split(url)[-1]

    if os.path.exists(dir_name):
        print('Directory already exists (%s)' % dir_name)
    else:
        print('Downloading (%s)' % filename)

        while not isfileok(filename, cache_dir):
            print('Fetching file %s' % filename)
            download(url, filename, cache_dir)

        print('Unpacking file %s' % filename)
        stat = unzip('%s/%s' % (cache_dir, filename), dir_name)
        if stat[0] != 0:
            print('Error unpacking %s' % filename)
            return False

        if post.has_key(dep):
            print('Running POST scripts')
            (status, out) = sh(post[dep])
            if status != 0:
                print out

        print('Successfully installed to %s' % dir_name)
        return True

if len(sys.argv) == 1:
    print_usage()
    sys.exit(1)

if match_command(sys.argv[-1]) == 'clean':
    print('Cleaning downloaded 3rd-party dependencies')
    sh('rm -r 3rd-party/ugm')
    sh('rm -r 3rd-party/maxflow')
    sh('rm -r 3rd-party/max-flow/maxflow-v3.0')
    sh('rm -r 3rd-party/dde-prnu')
    sys.exit(0)

if match_command(sys.argv[-1]) == 'dependencies':

    print('Configuring 3rd-party dependencies')

    if not os.path.exists('3rd-party'):
        print('Creating directory for 3rd party libraries (./3rd-party)')
        os.mkdir('3rd-party')

    for dep, url in deps.iteritems():
        filename = os.path.split(url)[-1]
        fetch_dependency(dep, url, '3rd-party')

    sys.exit(0)

if match_command(sys.argv[-1]) == 'data:maps':

    url = 'http://kt.agh.edu.pl/~korus/files/data-tifs-2016-maps.zip'
    fetch_dependency('tifs-2016-maps', url, 'data')

    sys.exit(0)

if match_command(sys.argv[-1]) == 'data:images':

    url = 'http://kt.agh.edu.pl/~korus/files/realistic-tampering-dataset.zip'
    fetch_dependency('images', url, 'data')

    sys.exit(0)

print('Unknown or ambiguous command: %s!' % sys.argv[-1])
sys.exit(2)
