#!/usr/bin/env python

import sys, os, errno, math
from collections import defaultdict, OrderedDict, Counter

def which(p):
    if os.path.isdir(p):
        return 'dir'
    if os.path.isfile(p):
        return 'file'

class _Seq:
    def __init__(self, tag, seq, quality=None):
        self.tag = tag
        self.seq = seq
        if quality:
            self.quality = quality

class Tag(OrderedDict):
    def make(self, tag):
        cells = tag.strip().split('|')
        cell_len = len(cells)
        n = 0
        while n + 1 < cell_len:
            self[cells[n]] = cells[n+1]
            n += 2

    def __str__(self):
        tags=[]
        for k in self:
            tags.extend([str(k), self[k]])
        return '|'.join(tags)

    def add(self, k, v):
        self[k] = v.replace('|','_')

    @classmethod
    def parse(cls, tag):
        ob = cls()
        ob.make(tag)
        return ob

class Seq(list):
    def _yield(self, ITER, infmt='fasta'):
        if infmt.lower() in ['fasta','fas','mfa','fna','fsa','fa']:
            tag, seq = '', ''
            for l in ITER:
                l = l.strip()
                if l and l[0] == '>':
                    if seq:
                        yield _Seq(tag.rstrip(), seq.strip().replace(' ',''))
                        seq = ''
                    tag = l[1:]
                else:
                    seq += l
            yield _Seq(tag.rstrip(), seq.strip().replace(' ',''))

    
    def _make(self, string, infmt='fasta'):
        if infmt.lower() in ['fasta','fas','mfa','fna','fsa','fa']:
            tag, seq = '', ''
            for l in string.splitlines():
                l = l.strip()
                if l and l[0] == '>':
                    if seq:
                        self.add(tag.rstrip(), seq.strip().replace(' ',''))
                        seq = ''
                    tag = l[1:]
                else:
                    seq += l
            self.add(tag.rstrip(), seq.strip().replace(' ',''))

        elif infmt.lower() == 'msf':
            _, alg = string.split('//')
            blocks = alg.strip().split('\n\n')
            tags = []
            tag2seq = defaultdict(str)
            for block in blocks:
                if not block:
                    continue
                for l in block.splitlines():
                    if not l:
                         continue
                    cs = l.split()
                    tag2seq[cs[0]] += ''.join(cs[1:])
                    if not cs[0] in tags:                
                        tags.append(cs[0])
            for tag in tags:
                self.add(tag, tag2seq[tag])

        elif infmt.lower() == 'phylip':
            i = 0
            for block in string.split('\n\n'):
                if not block:
                    continue
                j = 0
                for l in block.splitlines():
                    if not l:
                         continue
                    if l[0] != ' ':
                        cs = l.split()
                        self.add(cs[0], ''.join(cs[1:]))
                    elif i > 0:
                        self[j].seq += l.strip().replace(' ','')
                    j += 1
                i += 1

    def add(self, tag, seq, quality=None):
        self.append(_Seq(tag, seq, quality))
    
    def erase_common_gaps(self):
        pos2num  = defaultdict(int)
        for o in self:
            for idx in [i for i, x in enumerate(o.seq) if x == '-']:
                pos2num[idx] += 1
        ob_len = len(self)
        pos_to_erase = []
        for k, v in pos2num.iteritems():
            if v == ob_len:
                pos_to_erase.append(k)
        pos_to_erase.sort(reverse=True)
        for o in self:
            for pos in pos_to_erase:
                o.seq = o.seq[:pos] + o.seq[pos+1:]

    def upper(self):
        for o in self:
            if not o.seq.isupper():
                o.seq=o.seq.upper()

    def lower(self):
        for o in self:
            if not o.seq.islower():
                o.seq=o.seq.lower()

    def rm(self, char):
        for o in self:
            o.seq = o.seq.replace(char, '') 

    def type(self):
        seq_len = 0
        c = Counter(self[0].seq.lower())
        for k,v in c.iteritems():
            if k != '-':
                seq_len += v
        if (c['a'] + c['c'] + c['g'] + c['t'] + c['u']) / float(seq_len) < 0.4:
            return 'Protein'
        elif c['u']:
            return 'RNA'
        else:
            return 'DNA'

    def seq_len(self, raise_error=False):
        seqLen=0
        for ob in self:
            currentLen=len(ob.seq)
            if seqLen:
                if seqLen!=currentLen:
                    if raise_error:
                        sys.stderr.write('Mismatch in sequence len.\n')
                        sys.exit(0)
                    else:
                        return False
            else:
                seqLen=currentLen
        return seqLen

    def write(self, oh, outfmt='fasta', subformat='interleaved', max_tag_len=10, line_len=60, block_len=10, quiet=False):
        if outfmt.lower() in ['phylip']:
            seqs = []
            tags = []
            max_seq_len = 0
            seq_offset = max_tag_len + 3
            for o in self:
                tags.append(o.tag[:max_tag_len])
                curr_seq_len = len(o.seq)
                if curr_seq_len > max_seq_len:
                    max_seq_len = curr_seq_len
                count = 0
                seq = []
                while 1:
                    frag = o.seq[count:count+block_len]
                    if not frag:
                        break
                    seq.append(frag)
                    count += block_len
                seqs.append(seq)
            oh.write(' {} {}\n'.format(len(self), max_seq_len))
            block_num = int(line_len/block_len)
            for i in range(int(math.ceil(float(max_seq_len) / line_len))):
                for j in range(len(self)):
                    if i:
                        oh.write(' '*(seq_offset))
                    else:
                        oh.write(tags[j]+' '*(seq_offset-len(tags[j])))
                    oh.write(' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n')
                oh.write('\n')

        elif outfmt.lower() in ['nex', 'nxs', 'nexus']:
            seqs = []
            tags = []
            max_seq_len = 0
            seq_offset = max_tag_len + 3
            for o in self:
                tags.append(o.tag[:max_tag_len])
                curr_seq_len = len(o.seq)
                if curr_seq_len > max_seq_len:
                    max_seq_len = curr_seq_len
                count = 0
                seq = []
                while 1:
                    frag = o.seq[count:count+block_len]
                    if not frag:
                        break
                    seq.append(frag)
                    count += block_len
                seqs.append(seq)
            block_num = int(line_len/block_len)
            datatype=self.type().lower()
            if not datatype:
                sys.stderr.write('unknown datatype\n')
                sys.exit(0)
            oh.write('#NEXUS\n\nbegin data;\n')
            oh.write('\tdimensions ntax={} nchar={};\n'.format(len(self), max_seq_len))
            oh.write('\tformat datatype='+datatype+' interleave=yes gap=-;\n')
            oh.write('\tmatrix\n')
            for i in range(int(math.ceil(float(max_seq_len) / line_len))):
                for j in range(len(self)):
                    oh.write(tags[j]+' '*(13-len(tags[j]))+' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n')
                oh.write('\n')
            oh.write('\t;\nend;\n')

        elif outfmt.lower() in ['msf']:
            seqs = []
            tags = []
            max_seq_len = 0
            seq_offset = max_tag_len + 3
            for o in self:
                tags.append(o.tag[:tag_max_len])
                curr_seq_len = len(o.seq)
                if curr_seq_len > max_seq_len:
                    max_seq_len = curr_seq_len
                count = 0
                seq = []
                while 1:
                    frag = o.seq[count:count+block_len]
                    if not frag:
                        break
                    count += block_len
                    seq.append(frag.replace('-', '.'))
                seqs.append(seq)
            block_num = int(line_len/block_len)
            datatype=self.type()
            if datatype == 'Protein':
                datatype = 'P'
            else:
                datatype = 'N'

            oh.write('PileUp\n\n')
            oh.write(' MSF: {} TYPE: {} Check: 0 ..\n\n'.format(max_seq_len, datatype))
            for j in range(len(self)):
                oh.write(' Name: {} oo Len: {}  Check: 0 Weight: 10.00 ..\n'.format(tags[j], len(seqs[j])))
            oh.write('\n//\n\n')
            for i in range(int(math.ceil(float(max_seq_len) / line_len))):
                for j in range(len(self)):
                    oh.write(tags[j]+' '*(seq_offset-len(tags[j]))+' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n')
                oh.write('\n\n')

        elif outfmt.lower() in ['tsv']:
            for ob in self:
                oh.write(ob.tag+'\t'+ob.seq+'\n')

        elif outfmt.lower() in ['csv']:
            for ob in self:
                oh.write(ob.tag+','+ob.seq+'\n')

        elif outfmt.lower() in ['fasta','fas','mfa','fna']:
            for ob in self:                
                oh.write(">"+ob.tag+'\n')
                i = 0
                while 1:
                    frag = ob.seq[i:i+line_len]
                    if not frag:
                        break
                    oh.write(frag+'\n')
                    i += line_len
        
        if not quiet:
            sys.stdout.write('saved at {}\n'.format(oh.name))

    @classmethod
    def parse_string(cls, string, infmt='fasta'):
        o = cls()
        o._make(string, infmt=infmt)
        return o

    @classmethod
    def parse(cls, fh, infmt='fasta'):
        return cls.parse_string(fh.read(), infmt=infmt)

    @classmethod
    def iterparse(cls, fh, infmt='fasta'):
        o = cls()
        return o._yield(fh, infmt=infmt)

    @classmethod
    def parse_and_cath_files(cls, files):
        obs = []
        seq_len = None
        for p in files:
            s = Seq.parse(open(p))
            current_seq_len = len(s)
            if seq_len != None and current_seq_len != seq_len:
                sys.stderr.write('Mismatch in sequence # among files.\n')
                sys.exit(0)
            else:
                seq_len = current_seq_len
            obs.append(s)
        newseq = Seq()
        for i in range(seq_len):
            seq = ''
            for ob in obs:
                seq += ob[i].seq
            newseq.add(obs[0][i].tag, seq)
        return newseq

 
class Seq(Seq):
    def write_seq(self, write_ob, format, block_len, line_len, rm=None):
        if rm:
            self.rm(rm)
        #style seq output 
        if write_ob:
            if not format:
                format = args.o.name.split('.')[-1]
            self.write(write_ob, outfmt=format, block_len=block_len, line_len=line_len, quiet=False)
        else:
            self.write(sys.stdout, outfmt=format, block_len=block_len, line_len=line_len)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('a', metavar='action', choices=['cnvt', 'catv', 'cath'])
    parser.add_argument('i', metavar='in_file')
    parser.add_argument('-o', metavar='out_file', type=argparse.FileType('w'))
    parser.add_argument('-d', metavar='out_dir')
    parser.add_argument('-infmt', metavar='infmt', default='fasta')
    parser.add_argument('-inexts', metavar='inexts', nargs='+')
    parser.add_argument('-outfmt', metavar='outfmt', default='fasta')
    parser.add_argument('-ll', metavar='line_length', type=int, default=60)
    parser.add_argument('-bl', metavar='block_length', type=int, default=10)
    parser.add_argument('-rm')
    args = parser.parse_args()

    def mkdir_p(path):
        if path:
            try:    
                os.makedirs(path)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else:
                    raise

    class Files(list):
        def __init__(self, inpath, exts=[], depth=1):
            self.extend(self.make(inpath, exts = exts, depth = depth))

        def make(self, inpath, exts=[], depth = 1, depthCount = 1):
            if depthCount > depth:
                return None
            l=[]
            for filename in os.listdir(inpath):
                filepath = os.path.join(inpath, filename)
                if os.path.isdir(filepath):
                    l.extend(self.make(filepath, exts=exts, depth=depth, depthCount=depthCount+1))
                else:
                    if exts and not os.path.splitext(addr)[1] in exts:
                        continue
                    subpath = filename
                    stem = os.path.split(filepath)[0]
                    for c in range(depthCount-1):
                        stem, base = os.path.split(stem)
                        subpath = os.path.join(base, subpath)
                    l.append([filepath, subpath])
            return l

    if args.a == 'cnvt':
        #convert
        if which(args.i)=='dir':
            if not args.d:
                sys.stderr.write('A directory path for output files is required after \'-d\' when multiple input files are given.\n')
                sys.exit(0)
            mkdir_p(args.d)
            for f, subf in Files(args.i, exts=args.inexts, depth=float('inf')):
                stem, filename = os.path.split(subf)
                basename, ext = os.path.splitext(filename)
                ofh_dir = os.path.join(args.d, stem)
                mkdir_p(ofh_dir)
                ofh_name = os.path.join(ofh_dir, basename+'.'+args.outfmt)
                with open(ofh_name, 'w') as ofh:
                    Seq.parse(open(f), infmt=args.infmt).write_seq(ofh, args.outfmt, args.bl, args.ll, rm=args.rm)
        else:
            Seq.parse(open(args.i), infmt=args.infmt).write_seq(args.o, args.outfmt, args.bl, args.ll, rm=args.rm)
 
    if args.a == 'catv':
        if which(args.i)=='dir':
            s = Seq()
            for f, subf in Files(args.i, exts=args.inexts, depth=float('inf')):
                s.extend(Seq.parse(open(f)))
            s.write_seq(args.o, args.outfmt, args.bl, args.ll, rm=args.rm)
        
    if args.a == 'cath':
        if which(args.i)=='dir':
            seq = Seq.parse_and_cath_files([e[0] for e in Files(args.i, exts=args.inexts, depth=float('inf'))])
            seq.write_seq(args.o, args.outfmt, args.bl, args.ll, rm=args.rm)


