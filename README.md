# gmx_saltbridge

The gmx_saltbridge script allows for easy evaluation of internal salt bridge formation in a peptide or protein over a Gromacs molecular dynamics simulation trajectory. Charged residues are determined based on the primary sequence and the minimum distance between cationic and anionic groups are calculated for each frame of the trajectory. The salt bridge occupancy rate for each charge-pair is calculated based on how often the distance is less than a cutoff value and presented together with the mean distance in a csv file.     

**How to use**


The following dependencies are needed to correctly run the script: 

    Gromacs 

    Python3 with the NumPy and Pandas libraries 

  Running the script: 

    Download/copy the gmx_saltbridge.py script into the working directory for the simulation to be analysed 

    Create a simple text file containing only sequence of the peptide/protein in one letter code and name this file sequence.dat 

    Run the script in the terminal with the relevant arguments (seen using –h)  

    Results are shown in the terminal and saved to the results.txt file  


**Example output**

Following is an example output from the analysis of a solution simulation of the antimicrobial peptide Smp24. 

$ python3 gmx_saltbridge2.py -f traj.xtc -s md_3_500.tpr -b 0 -e 500000 -co 3.2 -ff Amber 

Output (occupancy in %, mean distance in nm): 

![image](https://user-images.githubusercontent.com/127429845/232488012-1e867cb1-7906-42fe-98ed-b260f391f7d9.png)

