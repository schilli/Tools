#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 

# additions:
# Copyright (c) 2005 Robert L. Campbell
from pymol import cmd,stored
import seq_convert

QuietException = parsing.QuietException
tmp_editor = "_tmp_editor"
tmp_ed_save = "_tmp_ed_save"
tpk1 = "_tmp_tpk1"
tpk2 = "_tmp_tpk2"

# attach_amino_acid is from pymol/modules/pymol/editor.py 
def attach_amino_acid(selection,amino_acid,phi,psi):
  if not selection in cmd.get_names("selections"):
    if amino_acid in cmd.get_names("objects"):
      print " Error: an object with than name already exists"
      raise QuietException
    cmd.fragment(amino_acid)
    if cmd.get_setting_legacy("auto_remove_hydrogens"):
      cmd.remove("(hydro and %s)"%amino_acid)
    if cmd.count_atoms("((%s) and name c)"%amino_acid,quiet=1):
      cmd.edit("((%s) and name c)"%amino_acid)
  else:
    cmd.fragment(amino_acid,tmp_editor)
    if cmd.count_atoms("((%s) and elem n)"%selection,quiet=1):
      cmd.select(tmp_ed_save,"(%s)"%selection)
      cmd.iterate("(%s)"%selection,"stored.resv=resv")
      stored.resi = str(stored.resv-1)
      cmd.alter(tmp_editor,"resi=stored.resi")
      cmd.fuse("(%s and name C)"%(tmp_editor),"(pk1)",2)
      if cmd.get_setting_legacy("auto_remove_hydrogens"):
        cmd.remove("(pkmol and hydro)")
      cmd.set_dihedral("(name ca and neighbor pk2)",
                            "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
      cmd.set_geometry("pk2",3,3) # make nitrogen planer
#      if ss:
      cmd.select(tpk1,"pk2")
      cmd.select(tpk2,"pk1")
      if amino_acid[0:3]!='pro':
        cmd.set_dihedral( # PHI
            "(name c and neighbor (name ca and neighbor "+tpk1+"))", # C
            "(name ca and neighbor "+tpk1+")", # CA 
            tpk1, # N
            tpk2, # C
            phi)
      cmd.set_dihedral( # PSI (n-1)
          tpk1, # N
          tpk2, # C
          "(name ca and neighbor "+tpk2+")", # CA
          "(name n and neighbor (name ca and neighbor "+tpk2+"))", # C
          psi)
      cmd.delete(tpk1)
      cmd.delete(tpk2)
      sele = ("(name N and (byres neighbor %s) and not (byres %s))"%
                (tmp_ed_save,tmp_ed_save))
      if cmd.count_atoms(sele,quiet=1):
        cmd.edit(sele)
      cmd.delete(tmp_ed_save)

    elif cmd.count_atoms("((%s) and elem c)"%selection,quiet=1):
      cmd.select(tmp_ed_save,"(%s)"%selection)
      cmd.iterate("(%s)"%selection,"stored.resv=resv")
      stored.resi = str(stored.resv+1)
      cmd.alter(tmp_editor,"resi=stored.resi")
      cmd.fuse("(%s and name N)"%(tmp_editor),"(pk1)",2)
      if cmd.get_setting_legacy("auto_remove_hydrogens"):
        cmd.remove("(pkmol and hydro)")
      cmd.set_dihedral("(name ca and neighbor pk2)",
                            "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
      cmd.set_geometry("pk1",3,3) # make nitrogen planer
#      if ss:
      cmd.select(tpk1,"pk1")
      cmd.select(tpk2,"pk2")
      if amino_acid[0:3]!='pro':
        cmd.set_dihedral( # PHI
            tpk2, # C
            tpk1, # N
            "(name ca and neighbor "+tpk1+")", # CA 
            "(name c and neighbor (name ca and neighbor "+tpk1+"))", # C
            phi)
      cmd.set_dihedral( # PSI (n-1)
          "(name n and neighbor (name ca and neighbor "+tpk2+"))", # C
          "(name ca and neighbor "+tpk2+")", # CA
          tpk2, # C
          tpk1, # N
          psi)
      cmd.delete(tpk1)
      cmd.delete(tpk2)
      sele = ("(name C and (byres neighbor %s) and not (byres %s))"%
                (tmp_ed_save,tmp_ed_save))
      if cmd.count_atoms(sele,quiet=1):
        cmd.edit(sele)
      cmd.delete(tmp_ed_save)
    elif cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
      print " Error: please pick a nitrogen or carbonyl carbon to grow from."
      cmd.delete(tmp_editor)
      raise QuietException

  cmd.delete(tmp_editor)


def build_seq_phi_psi(seq_phi_psi_filename,line_style=0):
  """
  usage: build_seq_phi_psi filename

  will build the above sequence from information in a file stored as:
  resname phi psi
  resname phi psi
  ...

  where resname is the 1-letter code for the amino acid

  The created object will be named for the first amino acid in the sequence,
  unless a pk1 selection exists, then it will build onto that atom.

  Default is now to automatically show the structure as sticks, call it with
  line_style=1 to revert to showing lines

  """

  line_style = int(line_style)
# read file of residue name with associated phi,psi values to build
  lines = open(seq_phi_psi_filename).readlines()
  sequence = ''
  phi_list = []
  psi_list = []
  sequence3 = []
  counter = 0
  for l in lines:
    if len(l.split()) == 3:
      counter += 1
      seq,phi,psi = l.split()
      if len(seq) == 3:
        sequence3.append(seq)
      elif len(seq) == 1:
        sequence += seq
      else:
        print " Error: problem with sequence file -- not 3-letter or 1-letter code for residue name"
        raise QuietException
    else:
      continue


    phi_list.append(float(phi))
    psi_list.append(float(psi))

  print "Number of amino acid/phi/psi values read: ",counter
  if counter == 0:
    print " Error: problem with sequence file -- there should be three columns: amino-acid phi psi"
    raise QuietException

  if len(sequence3) > 0:
    sequence = seq_convert.seq3_to_seq1(sequence3)
  print "Building sequence: ",sequence
  seq3=seq_convert.seq1_to_seq3(sequence)
  #print seq3[0].lower()

  if 'pk1' in cmd.get_names("selections"):
    obj='pk1'
  else:
    obj=seq[0:3]

  attach_amino_acid(obj,seq3[0].lower(),phi_list[0],psi_list[0])
  print seq3[0].lower(),phi_list[0],psi_list[0]
  for i in range(1,len(seq3)):
    aa = seq3[i].lower()
    attach_amino_acid('pk1',aa,phi_list[i],psi_list[i-1])
    print seq3[i].lower(),phi_list[i],psi_list[i-1]

# hide lines and show sticks (comment out the next two lines, if you don't like this).
  if line_style == 0:
    print "hiding lines for: ",seq3[0]
    print "showing sticks for: ",seq3[0]
    cmd.hide('lines',seq3[0])
    cmd.show('sticks',seq3[0])

cmd.extend('build_seq_phi_psi',build_seq_phi_psi)
