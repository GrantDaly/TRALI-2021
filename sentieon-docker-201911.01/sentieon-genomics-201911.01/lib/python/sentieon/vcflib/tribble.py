import collections
import io
import struct
import time

from .compat import *

__all__ = ['TribbleIndex']

class Header(object):
    __slots__ = ('magic', 'type', 'version', 'filename', 'filesize',
        'timestamp', 'md5', 'flags', 'properties')
    def __init__(self, *args):
        for k,v in zip(self.__slots__, args):
            setattr(self, k, v)

class LinearIndex(object):
    __slots__ = ('chrom', 'width', 'longest', 'count', 'blocks')
    def __init__(self, *args):
        for k,v in zip(self.__slots__, args):
            setattr(self, k, v)

class TribbleIndex(object):
    MAGIC = 0x58444954
    INDEX_TYPE_LINEAR = 1
    INDEX_TYPE_INTERVAL_TREE = 2
    VERSION = 3
    SEQUENCE_DICTIONARY_FLAG = 0x8000

    DEFAULT_INDEX_BIN_WIDTH = 8000
    GVCF_INDEX_BIN_WIDTH = 128000
    MAX_FEATURES_PER_BIN = 100
    MAX_FEATURES_PER_INTERVAL = 600

    def __init__(self, idxf, mode='r'):
        if not idxf.endswith('.idx'):
            raise ValueError('File name suffix is not .idx')
        self.path = idxf
        self.mode = mode
        if mode[0:1] == 'r':
            self.load()
        elif mode[0:1] == 'w':
            self.init()
        else:
            raise IOError('Mode ' + mode + ' not supported')

    def load(self):
        self.indices = collections.OrderedDict()

        with io.open(self.path, 'rb') as fp:
            data = fp.read()
            off = 0

            s = struct.Struct('<iii')
            magic, type, version = s.unpack_from(data, off)
            off += s.size
            if (magic != self.MAGIC or version != self.VERSION or
                type != self.INDEX_TYPE_LINEAR):
                raise RuntimeError('Bad magic/type/version')

            file,_ = data[off:].split(b'\0',1)
            off += len(file) + 1
            file = file.decode()

            s = struct.Struct('<QQ')
            size,time = s.unpack_from(data, off)
            off += s.size

            md5,_ = data[off:].split(b'\0',1)
            off += len(md5) + 1

            s = struct.Struct('<ii')
            flags,nprop = s.unpack_from(data, off)
            off += s.size
            if nprop > 0:
                t = data[off:].split(b'\0', 2*nprop)
                if len(t) != 2*nprop+1:
                    raise RuntimeError('Incorrect property count')
                data = t.pop(); off = 0
                t = [s.decode() for s in t]
                properties = zip(*[iter(t)]*2)
            else:
                properties = []

            self.header = Header(magic, type, version, file,
                size, time, md5, flags, properties)

            s = struct.Struct('<i')
            nchrs, = s.unpack_from(data, off)
            off += s.size
            for _ in xrange(nchrs):
                chrom,_ = data[off:].split(b'\0',1)
                off += len(chrom) + 1
                chrom = chrom.decode()
                if type ==  self.INDEX_TYPE_LINEAR:
                    s = struct.Struct('<iiiii')
                    width,bins,longest,_,count = s.unpack_from(data, off)
                    off += s.size;
                    s = struct.Struct('<'+'Q'*(bins+1))
                    blocks = s.unpack_from(data, off)
                    off += s.size
                    ci = LinearIndex(chrom, width, longest, count, blocks)
                else:
                    ci = None
                self.indices[chrom] = ci

    def save(self):
        if self.header is None:
            return
        self.add(None, 0, 0, self.end)
        h = self.header
        h.filesize = self.end
        h.timestamp = int(time.time())
        with io.open(self.path, 'wb') as fp:
            fp.write(struct.pack('<iii', h.magic, h.type, h.version))
            fp.write(h.filename.encode()); fp.write(b'\0')
            fp.write(struct.pack('<QQ', h.filesize, h.timestamp))
            fp.write(h.md5); fp.write(b'\0')
            fp.write(struct.pack('<ii', h.flags, len(h.properties)))
            for k,v in h.properties:
                fp.write(k.encode()); fp.write(b'\0')
                fp.write(v.encode()); fp.write(b'\0')
            fp.write(struct.pack('<i', len(self.indices)))
            for k,v in iteritems(self.indices):
                fp.write(k.encode()); fp.write(b'\0')
                fp.write(struct.pack('<iiiii', v.width,
                    len(v.blocks)-1, v.longest, 0, v.count))
                fp.write(struct.pack('<'+'Q'*len(v.blocks), *v.blocks))
        self.header = None

    def query(self, c, s, e):
        ci = self.indices.get(c)
        if ci is None:
            return []
        s = max(s, ci.longest) - ci.longest
        i = s // ci.width
        if i >= len(ci.blocks):
            return []
        return [(ci.blocks[i], ci.blocks[-1])]

    @staticmethod
    def merge(ranges, shift):
        p = None
        for r in sorted(ranges):
            if p is None:
                p = r
            elif r[0] >> shift > p[1] >> shift:
                yield p
                p = r
            else:
                p = (p[0], max(p[1],r[1]))
        if p is not None:
            yield p

    def init(self):
        self.header = Header(self.MAGIC, self.INDEX_TYPE_LINEAR, self.VERSION,
            self.path[:-4], 0, 0, b'', 0, [])
        self.width = (self.path.endswith('.g.vcf.idx')
            and self.GVCF_INDEX_BIN_WIDTH or self.DEFAULT_INDEX_BIN_WIDTH)
        self.indices = collections.OrderedDict()
        self.ci = None
        self.pos = 0
        self.end = 0

    def add(self, c, s, e, off):
        if self.ci and self.ci.chrom != c:
            self.ci.blocks.append(self.end)
            self.optimize(self.ci)
            self.ci = None
        if self.ci is None and c is not None:
            self.ci = LinearIndex(c, self.width, 0, 0, [])
            self.indices[c] = self.ci
            self.pos = 0
        if self.ci:
            assert self.ci.chrom == c and s >= self.pos
            bin = s // self.ci.width
            if bin >= len(self.ci.blocks):
                self.ci.blocks += [self.end] * (bin+1-len(self.ci.blocks))
            self.ci.longest = max(self.ci.longest, e-s)
            self.ci.count += 1
            self.pos = s
        self.end = off

    def optimize(self, ci):
        max_density = self.MAX_FEATURES_PER_BIN
        maxsize = max(ci.blocks[i] - ci.blocks[i-1]
            for i in xrange(1, len(ci.blocks)))
        fullsize = ci.blocks[-1] - ci.blocks[0]
        scale = (max_density * fullsize) // (ci.count * maxsize)
        if scale > 1:
            bins = (len(ci.blocks)-1 + scale-1) // scale
            ci.blocks = [ci.blocks[i*scale] for i in xrange(bins)]
            ci.width *= scale

# vim: ts=4 sw=4 expandtab
