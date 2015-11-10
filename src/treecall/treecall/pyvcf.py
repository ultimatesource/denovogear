#!/usr/bin/env python

# Author: Ni Huang <nihuang at genetics dot wustl dot edu>

from __future__ import print_function
import sys
from os.path import splitext
import gzip
from memoize import memoized

def warning(*obj):
    print(*obj, file=sys.stderr)


@memoized
def PL10_order(ref, alt):
    if alt == '':
        a = [ref]
    else:
        a = [ref] + alt.split(',')
    aidx = {a[i]:i for i in xrange(0,len(a))}
    nt = ('A','C','G','T')
    for x in nt:
        if x not in aidx:
            aidx[x] = len(a)
    return [aidx[nt[i]]+aidx[nt[j]]*(aidx[nt[j]]+1)/2 for i in xrange(4) for j in xrange(i,4)]


@memoized
def DPR4_order(ref, alt):
    if alt == '':
        a = [ref]
    else:
        a = [ref] + alt.split(',')
    aidx = {a[i]:i for i in xrange(0,len(a))}
    nt = ('A','C','G','T')
    for x in nt:
        if x not in aidx:
            aidx[x] = len(a)
    return [aidx[x] for x in nt]


class Vcf(object):
    def __init__(self, line=None, fixed_only=False):
        if line is not None:
            self.line = line.rstrip('\n')
            self.lazy_parse(fixed_only)

    def lazy_parse(self, fixed_only=True):
        fields = self.line.split('\t')
        self.CHROM = fields[0]
        self.POS = fields[1]
        self.ID = fields[2]
        self.REF = fields[3]
        self.ALT = fields[4].rstrip('<X>').rstrip(',')
        self.QUAL = fields[5]
        self.FILTER = fields[6]
        self.INFO = fields[7]
        self.info = None
        if not fixed_only:
            self.FMT = fields[8]
            self.extra = fields[9:]
            self.gtypes = {}

    def extract_gtype(self, tag, fmt, func=None, *args):
        if tag in self.gtypes:
            return self.gtypes[tag]
        gtype = []
        for s in self.extra:
            content = (s.split(':'))[fmt[tag]]
            if func is None:
                gtype.append(content)
            else:
                gtype.append(func(content, *args))
        self.gtypes[tag] = gtype
        return gtype

    def extract_info(self, tag=None, func=None, *args):
        if self.info is None:
            info = {}
            for s in self.INFO.split(';'):
                k,eq,v = s.partition('=')
                if eq == '':
                    v = True
                else:
                    v = v.split(',')
                    if len(v) == 1:
                        v = v[0]
                info[k] = v
            self.info = info
        if tag is None:
            return self.info
        else:
            value = self.info.get(tag)
            if func is not None:
                return func(value, *args)
            else:
                return value

    def get_PL10(self, x):
        order = PL10_order(self.REF, self.ALT)
        pl = x.split(',')
        n = len(pl)
        return [pl[i] if i < n else '255' for i in order]

    def get_DPR4(self, x):
        order = DPR4_order(self.REF, self.ALT)
        dpr = x.split(',')
        if max(order) >= len(dpr):
            warning(self)
        return [dpr[i] for i in order]

    def __str__(self):
        return self.line


class VcfFile(object):
    def __init__(self, filename=None, ifmt=None):
        if filename is not None:
            self.filename = filename
            self.ifmt = ifmt or self._get_format(filename)
            self.opened = False
            self.vcf_hdr = []
            self.samples = []
            self.fmt = {}
            self.read_vcf_header()
            if len(self.samples) > 0:
                self.read_FMT()

    def __iter__(self):
        self.open()
        if self.seekable:
            line = '#'
        else:
            line = self.fmt_line
        while True:
            if line[0] == '#':
                pass
            else:
                yield Vcf(line)
            try:
                line = self.next()
            except:
                self.close()
                break

    @staticmethod
    def parse_FMT(fmt):
        tags = fmt.split(':')
        return {v:i for i,v in enumerate(tags)}

    def _get_format(self, filename=None):
        filename = filename or self.filename
        base, ext = splitext(filename)
        return 'z' if ext == '.gz' else 'v'

    def open(self, force=False):
        if self.opened:
            if not force:
                return self.f
            else:
                warning('re-opening %s' % self.filename)
                self.close()
        self.seekable = True
        if self.ifmt == 'z':
            f = gzip.open(self.filename, 'rb')
        else:
            if self.filename == '-':
                f = sys.stdin
                self.seekable = False
            else:
                f = open(self.filename, 'r')
        self.opened = True
        self.f = f

    def close(self):
        if self.opened:
            self.f.close()
            self.opened = False

    def next(self):
        return self.f.next()

    def read_vcf_header(self):
        self.open(force=True)
        for line in self.f:
            if line[0] == '#':
                self.vcf_hdr.append(line.rstrip('\n'))
                if line[1] != '#':
                    fields = line.rstrip('\n').split('\t')
                    self.samples = fields[9:]
                    if self.seekable:
                        self.data_pos = self.f.tell()
                    break
            else:
                raise ValueError('No header present.')

    def skip_header(self):
        if not self.opened:
            raise ValueError('I/O operation on closed file.')
        if self.seekable:
            self.f.seek(self.data_pos)
        else:
            pass

    def read_FMT(self):
        if not self.opened:
            self.open(force=True)
        while True:
            line = self.next()
            if line[0] != '#':
                break
        self.fmt_line = line.rstrip('\n')
        fields = self.fmt_line.split('\t')
        self.fmt = VcfFile.parse_FMT(fields[8])
        if self.seekable:
            self.f.seek(self.data_pos)
