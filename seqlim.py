#!/usr/bin/env python

import sys
import os
import errno
import math
from collections import defaultdict, OrderedDict, Counter


class Path:
    @staticmethod
    def mkdir_p(path):
        if path:
            try:
                os.makedirs(path)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
                else:
                    raise

    @staticmethod
    def listfiles_r(inpath, exts=[], depth_count=1):
        for filename in os.listdir(inpath):
            filepath = os.path.join(inpath, filename)
            if os.path.isdir(filepath):
                for f in Path.listfiles_r(
                    filepath, exts=exts, depth_count=depth_count+1
                ):
                    yield f
            elif os.path.isfile(filepath):
                if exts and not os.path.splitext(filepath)[1] in exts:
                    continue
                subpath = filename
                stem = os.path.split(filepath)[0]
                for c in range(depth_count-1):
                    stem, base = os.path.split(stem)
                    subpath = os.path.join(base, subpath)
                yield filepath, subpath


class _Seq:
    def __init__(self, tag, seq, quality=None):
        self.tag = tag
        self.seq = seq


class Tag(OrderedDict):
    def make(self, tag):
        cells = tag.strip().split('|')
        cell_len = len(cells)
        n = 0
        while n + 1 < cell_len:
            self[cells[n]] = cells[n+1]
            n += 2

    def __str__(self):
        tags = []
        for k in self:
            tags.extend([str(k), self[k]])
        return '|'.join(tags)

    def add(self, k, v):
        self[k] = v.replace('|', '_')

    @classmethod
    def parse(cls, tag):
        ob = cls()
        ob.make(tag)
        return ob


class Seq(list):
    def parse_fasta(self, string):
        tag, seq = '', ''
        for l in string.splitlines():
            if not l:
                continue
            if l[0] == '>':
                if seq:
                    self.add(tag.strip(), seq.strip().replace(' ', ''))
                    seq = ''
                tag = l[1:]
            else:
                seq += l
        self.add(tag.strip(), seq.strip().replace(' ', ''))

    def parse_phylip(self, string):
        blocks = string.split('\n\n')
        i = 0
        for block in blocks:
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
                    self[j].seq += l.strip().replace(' ', '')
                j += 1
            i += 1

    def parse_msf(self, string):
        _, alg = string.split('//')
        blocks = alg.strip().split('\n\n')
        tags = []
        tag2seq = {}
        for block in blocks:
            if not block:
                continue
            for l in block.splitlines():
                if not l:
                    continue
                cs = l.split()
                if not cs[0] in tags:
                    tags.append(cs[0])
                    tag2seq[cs[0]] = ''
                tag2seq[cs[0]] += ''.join(cs[1:])
        for tag in tags:
            self.add(tag, tag2seq[tag])

    def make(self, string, infmt='fasta'):
        infmt = infmt.lower()
        if infmt in ('fasta', 'fas', 'mfa', 'fna', 'fsa', 'fa'):
            self.parse_fasta(string)
        elif infmt == 'msf':
            self.parse_msf(string)
        elif infmt == 'phylip':
            self.parse_phylip(string)
 
    def _yield(self, ITER, infmt='fasta'):
        if infmt.lower() in ('fasta', 'fas', 'mfa', 'fna', 'fsa', ' fa'):
            tag, seq = '', ''
            for l in ITER:
                l = l.strip()
                if l and l[0] == '>':
                    if seq:
                        yield _Seq(tag.rstrip(), seq.strip().replace(' ', ''))
                        seq = ''
                    tag = l[1:]
                else:
                    seq += l
            yield _Seq(tag.rstrip(), seq.strip().replace(' ', ''))

    def add(self, tag, seq):
        self.append(_Seq(tag, seq))
    
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
                o.seq = o.seq.upper()

    def lower(self):
        for o in self:
            if not o.seq.islower():
                o.seq = o.seq.lower()

    def rm(self, char):
        for o in self:
            o.seq = o.seq.replace(char, '') 

    def type(self):
        seq_len = 0
        c = Counter(self[0].seq.lower())
        for k, v in c.iteritems():
            if k != '-':
                seq_len += v
        if (c['a'] + c['c'] + c['g'] + c['t'] + c['u']) / float(seq_len) < 0.4:
            return 'Protein'
        elif c['u']:
            return 'RNA'
        else:
            return 'DNA'

    def seq_len(self, raise_error=False):
        seqLen = 0
        for ob in self:
            currentLen = len(ob.seq)
            if seqLen:
                if seqLen != currentLen:
                    if raise_error:
                        sys.stderr.write('Mismatch in sequence len.\n')
                        sys.exit(0)
                    else:
                        return False
            else:
                seqLen = currentLen
        return seqLen

    def write_phylip(
        self, oh,
        max_tag_len=10, line_len=60, block_len=10
    ):
        seq_len = len(self[0].seq)
        seq_num = len(self)
        seq_offset = max_tag_len + 3
        oh.write(' %s %s\n' % (seq_num, seq_len))
        FROM = 0
        for i in range(int(np.ceil(float(seq_len) / line_len))):
            for j in range(seq_num):
                if i:
                    oh.write(' '*seq_offset)
                else:
                    tag = self[j].tag[:max_tag_len]
                    oh.write(tag+' '*(seq_offset-len(tag)))
                seq = self[j].seq[FROM:FROM+line_len]
                oh.write(' '.join(self.chunks(seq, block_len))+'\n')
            FROM += line_len
            oh.write('\n')

    def write_fasta(self, oh, line_len=60):
        for ob in self:
            l = '>'+ob.tag+'\n'
            i = 0
            while 1:
                frag = ob.seq[i:i+line_len]
                if not frag:
                    break
                l += frag+'\n'
                i += line_len
            oh.write(l)


    def write(
        self, oh, outfmt='fasta',
        max_tag_len=10, line_len=60, block_len=10, quiet=True
    ):
        outfmt = outfmt.lower()
        if outfmt == 'phylip':
            tags, seqs = self.format_seqs(max_tag_len, line_len, block_len)
            seq_len = len(self[0].seq)
            seq_offset = max_tag_len + 3
            
            oh.write(' %s %s\n' % (len(self), seq_len))

            block_num = int(line_len/block_len)
            for i in range(int(math.ceil(float(seq_len) / line_len))):
                for j in range(len(self)):
                    if i:
                        oh.write(' '*seq_offset)
                    else:
                        oh.write(tags[j]+' '*(seq_offset-len(tags[j])))
                    oh.write(
                        ' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n'
                        )
                oh.write('\n')

        elif outfmt == 'tsv':
            for ob in self:
                oh.write(ob.tag+'\t'+ob.seq+'\n')

        elif outfmt == 'csv':
            for ob in self:
                oh.write(ob.tag+','+ob.seq+'\n')

        elif outfmt in ('fasta', 'fas', 'mfa', 'fna'):
            for ob in self:
                l = '>'+ob.tag+'\n'
                i = 0
                while 1:
                    frag = ob.seq[i:i+line_len]
                    if not frag:
                        break
                    l += frag+'\n'
                    i += line_len
                oh.write(l)

        elif outfmt in ('nex', 'nxs', 'nexus'):
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
            oh.write('\tdimensions ntax=%s nchar=%s;\n' % (len(self), max_seq_len))
            oh.write('\tformat datatype='+datatype+' interleave=yes gap=-;\n')
            oh.write('\tmatrix\n')
            for i in range(int(math.ceil(float(max_seq_len) / line_len))):
                for j in range(len(self)):
                    oh.write(tags[j]+' '*(13-len(tags[j]))+' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n')
                oh.write('\n')
            oh.write('\t;\nend;\n')

        elif outfmt == 'msf':
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
            oh.write(' MSF: %s TYPE: %s Check: 0 ..\n\n' % (max_seq_len, datatype))
            for j in range(len(self)):
                oh.write(
                    ' Name: %s oo Len: %s  Check: 0 Weight: 10.00 ..\n' % (tags[j], len(seqs[j]))
                )
            oh.write('\n//\n\n')
            for i in range(int(math.ceil(float(max_seq_len) / line_len))):
                for j in range(len(self)):
                    oh.write(tags[j]+' '*(seq_offset-len(tags[j]))+' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n')
                oh.write('\n\n')
        
        if not quiet:
            sys.stdout.write('saved at {}\n'.format(oh.name))

    def format_seqs(self, max_tag_len, line_len, block_len):
        tags, seqs = [], []
        for o in self:
            tags.append(o.tag[:max_tag_len])
            curr_seq_len = len(o.seq)
            if curr_seq_len > line_len:
                row_len = curr_seq_len
            count = 0
            chunks = []
            while 1:
                frag = o.seq[count:count+block_len]
                if not frag:
                    break
                chunks.append(frag)
                count += block_len
            tags.append(o.tag[:max_tag_len])
            seqs.append(chunks)
        return tags, seqs

    @classmethod
    def parse_string(cls, string, infmt='fasta'):
        o = cls()
        o.make(string, infmt=infmt)
        return o

    @classmethod
    def parse(cls, fh, infmt='fasta'):
        return cls.parse_string(fh.read(), infmt=infmt)

    @classmethod
    def iterparse(cls, fh, infmt='fasta'):
        o = cls()
        return o._yield(fh, infmt=infmt)

    @classmethod
    def cath_files(cls, files):
        obs = []
        seq_len = None
        for p in files:
            s = Seq.parse(open(p))
            current_seq_len = len(s)
            if seq_len is not None and current_seq_len != seq_len:
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
    def write_seq(self, outfile, outfmt, block_len, line_len, rm=None):
        if rm:
            self.rm(rm) # remove residues

        # style seq output 
        if outfile:
            if not outfmt:
                outfmt = args.o.name.split('.')[-1]
            with open(outfile, 'w') as ofh:
                self.write(
                    ofh, outfmt=outfmt,
                    block_len=block_len, line_len=line_len,
                    quiet=False
            )   
        else:
            self.write(
                sys.stdout, outfmt=outfmt,
                block_len=block_len, line_len=line_len
            )


if __name__ == '__main__':
    import argparse
    par = argparse.ArgumentParser()
    par.add_argument(
        'a', metavar='cnvt|catv|cath', 
        choices=['cnvt', 'catv', 'cath'],
        help="choose among 'cnvt' (convert), 'catv' (concatenate vertically), 'cath' (concatenate horizontally))"
    )
    par.add_argument(
        'i', metavar='PATH',
        help='MSA file or dir containing MSA files.'
    )
    par.add_argument(
        '-o', metavar='PATH',
        help='file or dir for saving MSA(s).'
    )
    par.add_argument(
        '-infmt', metavar='fasta|phylip|msf', default='fasta',
        help="input format. default=fasta"
    )
    par.add_argument(
        '-outfmt', metavar='fasta|phylip|nex|msf|tsv|csv', default='fasta',
        help="output format. default=fasta"
    )
    par.add_argument(
        '-line_length', metavar='INTEGER', type=int, default=60,
        help="sequence length of each line. default=60"
    )
    par.add_argument(
        '-block_length', metavar='INTEGER', type=int, default=10,
        help="the length of each block in lines. default=10"
    )
    par.add_argument(
        '-remove', metavar="S1",
        help="the specified residue will be removed."
    )
    par.add_argument('-inexts', metavar='inexts', nargs='+')
    args = par.parse_args()

    if args.a == 'cnvt': # convert
        if os.path.isdir(args.i):
            if not args.o:
                sys.stderr.write("a dir path for output files is needed.\n")
                exit()
            Path.mkdir_p(args.o)
            for f, subf in Path.listfiles_r(args.i, exts=args.inexts):
                stem, filename = os.path.split(subf)
                basename, ext = os.path.splitext(filename)
                ofh_dir = os.path.join(args.o, stem)
                Path.mkdir_p(ofh_dir)
                ofh_name = os.path.join(ofh_dir, basename+'.'+args.outfmt)
                
                with open(f) as ifh:
                    seq = Seq.parse(ifh, infmt=args.infmt)
                    seq.write_seq(
                        ofh_name, args.outfmt, args.block_length, args.line_length, rm=args.remove
                    )
        else:
            with open(args.i) as ifh:
                seq = Seq.parse(ifh, infmt=args.infmt)
                seq.write_seq(args.o, args.outfmt, args.block_length, args.line_length, rm=args.remove)

    elif args.a == 'catv' and os.path.isdir(args.i):
        s = Seq()
        for path, subpath in Path.listfiles_r(args.i, exts=args.inexts):
            with open(path) as f:
                s.extend(Seq.parse(f))
            s.write_seq(args.o, args.outfmt, args.block_length, args.line_length, rm=args.remove)
        
    elif args.a == 'cath' and os.path.isdir(args.i):
        seq = Seq.cath_files([e[0] for e in Path.listfiles_r(args.i, exts=args.inexts)])
        seq.write_seq(args.o, args.outfmt, args.block_length, args.line_length, rm=args.remove)

