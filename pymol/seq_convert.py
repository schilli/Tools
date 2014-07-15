#! /usr/bin/env python
# Copyright (c) 2004 Robert L. Campbell

import sys,re,getopt

#import MyPDB

def usage():
  print "usage: seq_convert.py <options> input_file"
  print "       - input_file may be a file name or left blank (therefore input comes from stdin)"
  print "       - output is written to standard out\n"
  print "       Options can be:"
  print "       -i, --in [1,3,pdb,pdbcoord]    input format: 1 or 3-letter code, or PDB file <default 1>"
  print "                            'pdbcoord' means force the use of the coordinates to extract the sequence"
  print "       -o, --out [1,3]       output format: 1 or 3-letter code <default 1>"
  print "       -r, --res [#]         where # is number of residues/line to be printed <default 80>\n"
  print "       --num         print sequence numbers at the end of each line <default>"
  print "       --nonum       do not print sequence numbers at the end of each line"
  print "       --nohead      do not print sequence file header"
  print "       --help        print this help message>"
  print "       -p, --pir     add PIR formatting"
  print "       -q, --quiet     don't write messages only sequence"
  sys.exit(0)

def get_options():
# defaults
  in_format = '1'
  out_format = '1'
  res_per_line = 80
  print_seq_num = 0
  print_seq_head = 1
  file_in = sys.stdin
  pir = 0
  quiet=0

  try:
    opts,args = getopt.getopt(sys.argv[1:],'f:hi:o:r:pq',['help','file=','in=','out=','res=','num','nonum','head','nohead','pir','quiet'])
  except:
    sys.stderr.write("\n*********************************\n")
    sys.stderr.write("\n      Unknown options %s\n" % str(sys.argv[1:]))
    sys.stderr.write("\n*********************************\n\n")
    usage()

  if len(args) == 1:
    file_in = args[0]

# not sure why this was here?
#  if len(opts) == 0:
#    usage()
#  else:
  for o,a in opts:
    if o in ('-h', '--help'):
      usage()
    elif o in ('-f','--file'):
      file_in = a
    elif o in ('-i', '--in'):
      in_format = a
    elif o in ('-o', '--out'):
      out_format = a.lower()
      if out_format not in ('1','3'):
        sys.stderr.write('\n\nERROR:\nOutput format specification not understood: %s\n\n' % out_format)
        sys.exit(1)
    elif o in ('-r', '--res'):
      res_per_line = int(a)
    elif o in ('--num'):
      print_seq_num = 1
    elif o in ('--nonum'):
      print_seq_num = 0
    elif o in ('--head'):
      print_seq_head = 1
    elif o in ( '--nohead'):
      print_seq_head = 0
    elif o in ('-p', '--pir'):
      pir = 1
    elif o in ('-q', '--quiet'):
      quiet = 1

  return file_in,in_format,out_format,res_per_line,print_seq_num,print_seq_head,pir,quiet

def seq1_to_seq3(seq1):
  """
  Convert array of 1-letter code sequence to 3-letter code.
  seq1 can be an array or a string.
  Return seq3 as array.
  """
  res3 = ['---','ala','asn','asp','arg','cys','gln','glu',
           'gly','his','ile','leu','lys','met','pro','phe','ser',
           'thr','trp','tyr','val','unk',
           'ALA','ASN','ASP','ARG','CYS','GLN','GLU',
           'GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER',
           'THR','TRP','TYR','VAL','UNK']
  res1 = '-andrcqeghilkmpfstwyvxANDRCQEGHILKMPFSTWYVX'

  if len(seq1) > 1:
    seq3 = []
    for a1 in seq1:
      try:
        a3 = res3[res1.index(a1)]
        seq3.append(a3)
      except ValueError, err:
        print "%s:  No match for residue: %s" % (err,a1)
  else:
    seq3 = res3[res1.index(seq1)]

  return seq3

def to_upper(a):
  lower='abcdefghijklmnopqrstuvwxyz'
  upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  b = ''
#  if type(a) == type('STRING'):
  for i in range(len(a)):
    try:
      b += upper[lower.index(a[i])]
    except:
      b += a[i]
  return b
#  else:
#    return a

def seq3_to_seq1(seq3):
  """
  Convert array of 3-letter code sequence to 1-letter code
  Return seq1 as string
  """
  seq1 = ''
  res3 = ['---','ala','asn','asp','arg','cys','gln','glu',
           'gly','his','ile','leu','lys','met','pro','phe','ser',
           'thr','trp','tyr','val','unk', 
           'ALA','ASN','ASP','ARG','CYS','GLN','GLU',
           'GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER',
           'THR','TRP','TYR','VAL','UNK',]
# wanted to include '^M' and '^L' in the res1 string to match 'XXX','YYY' below
#           'THR','TRP','TYR','VAL','UNK','XXX','YYY']
  res1 = '-andrcqeghilkmpfstwyvxANDRCQEGHILKMPFSTWYVX'

#force seq3 to be a list, in case it is just a single amino acid in 3-letter code
  if type(seq3) != type(list()):
    seq3 = [seq3]
  for a3 in seq3:
#    a3 = to_upper(a3)
    a3 = a3.upper()
    # strip trailing spaces in case the three letter code came from a MyPDB Protein sequence listing
    try:
      a1 = res1[res3.index(a3.strip())]
    except ValueError:
      a1 = 'x'
    seq1 += a1

  return seq1

def write_out(seq,blank):
  """
  seq is a list of residues (either one or 3 letter code)
  blank is a blank string to insert at the end of incomplete lines 
  (' ' for 1-letter code and '    ' for 3-letter code)
  """

  if len(seq)%res_per_line == 0:
    num_lines = len(seq)/res_per_line
  else:
    num_lines = len(seq)/res_per_line + 1

  for i in range(num_lines):
    out = ''
    j = i*res_per_line
    k = (i+1)*res_per_line

  # check for last line (might be incomplete)
    if k <= len(seq):              # not last line of sequence
      for aa in seq[j:k]:
        if out_format == '1':
          out = out + aa
        else:
          out = out + '%-4s' % aa
      # write sequence number of last residue in the line
      if print_seq_num:
        out = out + '   %5d' % k

    else:                          #last line of sequence
      for aa in seq[j:len(seq)]:
        if out_format == '1':
          out = out + aa
        else:
          out = out + '%-4s' % aa
        # now add blanks to fill in space
        num_blanks = res_per_line-(len(seq)-j)
      for n in range(num_blanks):
        out = out + blank
      # write sequence number of last residue in the line
      if print_seq_num:
        out = out + '   %5d' % len(seq)
    print out

def readlines(file_in):
  if file_in == sys.stdin:
    lines = sys.stdin.readlines()
  else:
    lines = open(file_in).readlines()

  return lines

def readseq(file_in,in_format,quiet):
  header = []

  if in_format == '1':
    lines = readlines(file_in)
    seq_in = ''
    for line in lines:
      # ignore leading lines of standard sequence format
      if line[0] != '>':

        # strip off sequence numbers, spaces and end-of-line characters
        line = re.sub('[0-9\s]','',line)

        # append residues to sequence list
# not sure why I thought this was necessary
        #for i in range(len(line)):
        #  seq_in += line[i]
        seq_in += line

      else:
        header.append(line[:-1])

  elif in_format == '3':
    lines = readlines(file_in)
    seq_in = []
    for line in lines:
      # ignore leading lines of standard sequence format
      if line[0] != '>':
        # strip off sequence numbers, but leave spaces
        line = re.sub('[0-9]','',line)
        # append residues to sequence list
        for res in line.split():
          seq_in.append(res)

      else:
        header.append(line[:-1])

  elif in_format == 'pdb' or in_format == 'pdbcoord':
    lines = readlines(file_in)
    seq_in = {}
    found = 0
    if in_format == 'pdb':
      # need to possibly separate chains.
      for line in lines:
        if line[0:6] == 'SEQRES':
          found = 1
          chain = line[11]
          words = line[19:71].split()
          for res in words:
            if chain in seq_in:
              seq_in[chain].append(res)
            else:
              seq_in[chain] = [res]

    if not found or in_format == 'pdbcoord':
      if not quiet:
        sys.stderr.write("Getting sequence from coordinates...\n")
      import MyPDB
      p = MyPDB.Protein()
      p.readPDBlines(lines)
      seq_in = p.sequence

  return header,seq_in

if __name__ == '__main__':
  file_in,in_format,out_format,res_per_line,print_seq_num,print_seq_head,pir,quiet = get_options()

  seq_out = []

  header,seq_in = readseq(file_in,in_format,quiet)

# print header back out
  if print_seq_head:
    for head in header:
      print head

  if pir and file_in != sys.stdin:
    print ">%s;" % file_in
# do conversion and write out
  if in_format == '1' and out_format == '3':
    seq_out = seq1_to_seq3(seq_in)
    write_out(seq_out,'    ')
  elif in_format == '3'  and out_format == '1':
    seq_out = seq3_to_seq1(seq_in)
    write_out(seq_out,' ')
  elif (in_format == 'pdb' or in_format == 'pdbcoord') and out_format == '3':
    for chain in seq_in.keys():
      seq_out = seq_in[chain]
      write_out(seq_out,'    ')
      print ""
  elif (in_format == 'pdb' or in_format == 'pdbcoord') and out_format == '1':
    for chain in seq_in.keys():
      seq_out = seq3_to_seq1(seq_in[chain])
      write_out(seq_out,' ')
      print ""
#  elif in_format == '1' and out_format == '1':
#    seq_out = seq_in
#    write_out(seq_out,' ')
#    print ""
#  elif in_format == '3' and out_format == '3':
#    seq_out = seq_in
#    write_out(seq_out,' ')
#    print ""
  else:
    seq_out = seq_in
    write_out(seq_out,' ')
    print ""

