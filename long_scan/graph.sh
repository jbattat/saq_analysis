rm list.txt
for filename in *.root
 
do python ~/Research/qpix-digital/prototype-software/analyze_root.py $filename
   echo $filename >> list.txt
done


   
		
	
