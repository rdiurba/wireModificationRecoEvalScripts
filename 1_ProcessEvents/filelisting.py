# Script taken from protoduneana used by ProtoDUNE (Developed by Tingjun Yang)
import os
from argparse import ArgumentParser as ap
# reused from protoduneana/michelremoving/Xcalo/filelisting.py
parser = ap()

parser.add_argument( "-r", type=str, help='Run number', default='4b')
#parser.add_argument( "-o", type=str, help='Time (s) to sleep between calls', default=)

args = parser.parse_args()

fout = open('file_list_' + args.r + '.txt', 'w')


#getnames = os.popen("samweb list-definition-files run4b_bnb_beam_off_pandora_reco2_reco2_reco2_all")
getnames = os.popen("samweb list-definition-files "+args.r)
filenames = getnames.readlines()
out_names = []
for fn in filenames:
  #print(fn)
  fn = fn.rstrip()
  fileloc = os.popen(f"samweb locate-file {fn}")
  fileloc = fileloc.readlines()
  fileloc = fileloc[0].split(":")
  fileloc = fileloc[1].split("(")
  
  #print (fileloc[0]+"/"+fn)
  print(fileloc[0].strip())
  out_names.append(fileloc[0].strip()+"/"+fn +"\n")

fout.writelines(out_names)
