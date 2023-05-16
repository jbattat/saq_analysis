import sys
import subprocess
import saq

flist = sys.argv[1]   # name of files containing a list of root files to analyze

files = saq.files_from_list(flist)
print(files)

for ff in files:
    subprocess.call(['python', 'analyze_root.py', ff])#, shell=True)

   
		
	
