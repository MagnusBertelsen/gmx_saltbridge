from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import subprocess
import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'



### Parse command line arguments ###
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

parser.add_argument("-f", "--trajectory_name", default="traj.xtc", help="Name of the simulation trajectory")
parser.add_argument("-s", "--tpr_name", default='topol.tpr', help="Name of production simluation tpr-file")
parser.add_argument("-b", "--time_start", default=0, type=int, help="Time for start of analysis in ps / time for once the peptide has reached full insertion into the bilayer in ps")
parser.add_argument("-e", "--end_time", default=500000, type=int, help="Time for end of analysis / end of simultion in ps")
parser.add_argument("-co", "--cutoff", default=0.32, type=float, help="Cutoff value for the distance at which a saltbridge is considered formed")
parser.add_argument("-ff", "--forcefield", default="none", help="Name of force field used in the simulation. If Amber or CHARMM is chosen, distances will be calculated based on the charged groups (anionic  = carboxylic acids (C-terminal, Aspartic acid, Glutamic acid), cationic = amine (N-terminal & Lysine) and guanidino groups (Arginine)). If none is selected, distances are calculated based on sidechains as a whole")


args = vars(parser.parse_args())



### Set up parameters ###

trajectory_name=args["trajectory_name"]
tpr_name=args["tpr_name"]
time_start=args["time_start"]
end_time=args["end_time"]
cutoff=args["cutoff"]
ff=args["forcefield"]



### Functions ###


def get_sequence():
  with open('sequence.dat', 'r') as file:
      sequence_raw = file.read().rstrip()

  sequence = []

  for res in sequence_raw:
      sequence.append(res)

  sequence = pd.DataFrame(sequence, columns=['Residue'])
  return sequence



def find_cat():

  sequence = get_sequence()
  
  res_nr = []
  res_type = []
  
  res_nr.append(1)
  res_type.append('N-terminal')
  
  for i in range(0,len(sequence)):
    if sequence.iat[i,0] == 'K' or sequence.iat[i,0] == 'R':
      res_nr.append(i+1)
      res_type.append(sequence.iloc[i,0])

  res_cat = pd.DataFrame(data = zip(res_nr,res_type), columns = ['Residue number', 'Residue type'])
  
  return res_cat
  

def find_an():

  sequence = get_sequence()
  
  res_nr = []
  res_type = []
  
  
  for i in range(0,len(sequence)):
    if sequence.iat[i,0] == 'D' or sequence.iat[i,0] == 'E':
      res_nr.append(i+1)
      res_type.append(sequence.iloc[i,0])
      
  res_nr.append(len(sequence))
  res_type.append('C-terminal')
  
  res_an = pd.DataFrame(data = zip(res_nr,res_type), columns = ['Residue number', 'Residue type'])
  
  return res_an


def gen_pairdist_input():

  res_cat = find_cat()
  res_an = find_an()
  
  f= open("input_index.txt","w+")
  f.write("del10-50"+'\n')
  
  if ff == "Amber": #Amber 
    if len(res_cat) > len(res_an):
     
     for i in range(0,len(res_an)):
      
      if res_an.iat[i,1] == 'C-terminal':
        f.write(f"a name OC1 | a name OC2 & r{res_an.iat[i,0]}"+'\n')     
      elif res_an.iat[i,1] == 'E':
        f.write(f"a name OE1 | a name OE2 & r{res_an.iat[i,0]}"+'\n')
      else:
        f.write(f"a name OD1 | a name OD2 & r{res_an.iat[i,0]}"+'\n')
        
     for i in range(0,len(res_cat)):
      if res_cat.iat[i,1] == 'N-terminal':
        f.write(f"a name N | a name H1 | a name H2 | a name H3 & r{res_cat.iat[i,0]}"+'\n')
      elif res_cat.iat[i,1] == 'K':
        f.write(f"a name NZ | a name HZ1 | a name HZ2 | a name HZ3 & r{res_cat.iat[i,0]}"+'\n')
      else:
        f.write(f"a name NH1 | a name NH2 | a name HH11 | a name HH12 | a name HH21 | a name HH22 & r{res_cat.iat[i,0]}"+'\n')
        
             
    else:
    
     for i in range(0,len(res_cat)):
      if res_cat.iat[i,1] == 'N-terminal':
        f.write(f"a name N | a name H1 | a name H2 | a name H3 & r{res_cat.iat[i,0]}"+'\n')
      elif res_cat.iat[i,1] == 'K':
        f.write(f"a name NZ | a name HZ1 | a name HZ2 | a name HZ3 & r{res_cat.iat[i,0]}"+'\n')
      else:
        f.write(f"a name NH1 | a name NH2 | a name HH11 | a name HH12 | a name HH21 | a name HH22 & r{res_cat.iat[i,0]}"+'\n')
    
     for i in range(0,len(res_an)):
      if res_an.iat[i,1] == 'C-terminal':
        f.write(f"a name OC1 | a name OC2 & r{res_an.iat[i,0]}"+'\n')     
      elif res_an.iat[i,1] == 'E':
        f.write(f"a name OE1 | a name OE2 & r{res_an.iat[i,0]}"+'\n')
      else:
        f.write(f"a name OD1 | a name OD2 & r{res_an.iat[i,0]}"+'\n')
  
  
  elif ff == "CHARMM": #CHARMM
    if len(res_cat) > len(res_an):
     
     for i in range(0,len(res_an)):
      
      if res_an.iat[i,1] == 'C-terminal':
        f.write(f"a name OT1 | a name OT2 & r{res_an.iat[i,0]}"+'\n')     
      elif res_an.iat[i,1] == 'E':
        f.write(f"a name OE1 | a name OE2 & r{res_an.iat[i,0]}"+'\n')
      else:
        f.write(f"a name OD1 | a name OD2 & r{res_an.iat[i,0]}"+'\n')
        
     for i in range(0,len(res_cat)):
      if res_cat.iat[i,1] == 'N-terminal':
        f.write(f"a name N | a name H1 | a name H2 | a name H3 & r{res_cat.iat[i,0]}"+'\n')
      elif res_cat.iat[i,1] == 'K':
        f.write(f"a name NZ | a name HZ1 | a name HZ2 | a name HZ3 & r{res_cat.iat[i,0]}"+'\n')
      else:
        f.write(f"a name NH1 | a name NH2 | a name NH11 | a name NH12 | a name NH21 | a name NH22 & r{res_cat.iat[i,0]}"+'\n')
        
             
    else:
    
     for i in range(0,len(res_cat)):
      if res_cat.iat[i,1] == 'N-terminal':
        f.write(f"a name N | a name H1 | a name H2 | a name H3 & r{res_cat.iat[i,0]}"+'\n')
      elif res_cat.iat[i,1] == 'K':
        f.write(f"a name NZ | a name HZ1 | a name HZ2 | a name HZ3 & r{res_cat.iat[i,0]}"+'\n')
      else:
        f.write(f"a name NH1 | a name NH2 | a name NH11 | a name NH12 | a name NH21 | a name NH22 & r{res_cat.iat[i,0]}"+'\n')
    
     for i in range(0,len(res_an)):
      if res_an.iat[i,1] == 'C-terminal':
        f.write(f"a name OT1 | a name OT2 & r{res_an.iat[i,0]}"+'\n')     
      elif res_an.iat[i,1] == 'E':
        f.write(f"a name OE1 | a name OE2 & r{res_an.iat[i,0]}"+'\n')
      else:
        f.write(f"a name OD1 | a name OD2 & r{res_an.iat[i,0]}"+'\n') 
  
  
  else: # General / no force field   
   if len(res_cat) > len(res_an):
     
    for i in range(0,len(res_an)):
      if res_an.iat[i,1] == 'C-terminal':
        f.write(f"5 & r{res_an.iat[i,0]}"+'\n')
      else:
        f.write(f"8 & r{res_an.iat[i,0]}"+'\n')
        
    for i in range(0,len(res_cat)):
      if res_cat.iat[i,1] == 'N-terminal':
        f.write(f"5 & r{res_cat.iat[i,0]}"+'\n')
      else:
        f.write(f"8 & r{res_cat.iat[i,0]}"+'\n')
        
             
   else:
    
    for i in range(0,len(res_cat)):
      if res_cat.iat[i,1] == 'N-terminal':
        f.write(f"5 & r{res_cat.iat[i,0]}"+'\n')
      else:
        f.write(f"8 & r{res_cat.iat[i,0]}"+'\n')
    
    for i in range(0,len(res_an)):
      if res_an.iat[i,1] == 'C-terminal':
        f.write(f"5 & r{res_an.iat[i,0]}"+'\n')
      else:
        f.write(f"8 & r{res_an.iat[i,0]}"+'\n')
       
          
  f.write("del0-9"+'\n')
  f.write("q"+'\n')
    
  f.close()
  
  
  
  if len(res_cat) > len(res_an):  
      
    for i in range(0, len(res_an)):
      f= open(f"input{i}.txt","w+")
      
      f.write(f"{i}"+'\n')
      
      for x in range(0, len(res_cat)):
        f.write(f"{x+len(res_an)}"+'\n')  
      
      f.close()
      
  else:
    
    for i in range(0, len(res_cat)):
      f= open(f"input{i}.txt","w+")
      
      f.write(f"{i}"+'\n')
      
      for x in range(0, len(res_an)):
        f.write(f"{x+len(res_an)}"+'\n')  
      
      f.close()        
  
  
def pairdist(trajectory_name, tpr_name, time_start, end_time, cutoff):
  
  res_cat = find_cat()
  res_an = find_an()  
    
  f= open("pairdist.sh","w+")
  f.write(
  "#!/usr/bin/bash"+'\n'
  
  "tpr_name=${1}"+'\n'
  "xtc_name=${2}"+'\n'
  "start_time=${3}"+'\n'
  "end_time=${4}"+'\n'
    
  "gmx make_ndx -f ${1} -o pairdist.ndx < input_index.txt"+'\n')
    
  if len(res_cat) > len(res_an):     
  
    for i in range(0, len(res_an)):
      f.write(f"gmx pairdist -f ${2} -s ${1} -n pairdist.ndx -type min -o pairdist{i}.xvg -xvg none -b ${3} -e ${4} -cutoff 0.32 < input{i}.txt"+'\n')

  f.close()  

  os.chmod("pairdist.sh", 0o0777)
  
  subprocess.run(["./pairdist.sh", str(tpr_name), str(trajectory_name), str(time_start), str(end_time)])
  
  os.remove("pairdist.sh")
    
    

def analyse_single_residue(short_ion_df, long_ion_df, pairdist_file_name, ion_number, cutoff):
  
  names = []
  
  names.append('Time')
  
  ion = str(short_ion_df.iat[ion_number,1]) + str(short_ion_df.iat[ion_number,0])
  	
  for i in range (0, len(long_ion_df)):
    x = f"{ion}_{str(long_ion_df.iat[i,1]) + str(long_ion_df.iat[i,0])}"
    names.append(x)
  
  data = pd.read_csv(pairdist_file_name, delim_whitespace=True, header=None, names = names)
  
  occupancy = []
  
  for i in range(1, len(names)):
    prop = data.iloc[:, i].value_counts(normalize = True)
    occupancy.append((1-prop[0.32])*100)
  
  mean = data[data < cutoff]
  mean = mean.mean()
  mean = mean.drop(labels = 'Time')
  
  names.pop(0)
  
  results = pd.DataFrame(data = [names, occupancy, mean])
  results = results.T
  results = results.rename(columns = {0:'Residues', 1:'Occupancy rate', 2:'Mean distance'}) 
  
  return results
  

def analysis_full():
  
  res_an = find_an()
  res_cat = find_cat()
  
  if len(res_an) < len(res_cat):
  
    results = pd.DataFrame(columns = ['Residues', 'Occupancy rate', 'Mean distance'])
    
    for i in range(0, len(res_an)):
      pairdist_file_name = f"pairdist{i}.xvg"
      results = pd.concat([results, analyse_single_residue(res_an, res_cat, pairdist_file_name, i, cutoff)])
  
    results.to_csv('results.txt', index = False)
  
  else:
  
    results = pd.DataFrame(columns = ['Residues', 'Occupancy rate', 'Mean distance'])
    
    for i in range(0, len(res_cat)):
      pairdist_file_name = f"pairdist{i}.xvg"
      results = pd.concat([results, analyse_single_residue(res_cat, res_an, pairdist_file_name, i, cutoff)])
  
    results.to_csv('results.txt', index = False)
  
  return results
  






def main(trajectory_name, tpr_name, time_start, end_time, cutoff):
  

  gen_pairdist_input()
  pairdist(trajectory_name, tpr_name, time_start, end_time, cutoff)
  analysis_full()
  print(pd.read_csv('results.txt'))




if __name__ == '__main__':
  main(trajectory_name, tpr_name, time_start, end_time, cutoff)
  
  
  
