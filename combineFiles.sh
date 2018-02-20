#merges a tracking file (INFILE) with another combined file (OUTFILE)
#if outfile does not exist, it is created

INFILE=$1
OUTFILE=$2

#Get number of fields
   N_FIELDS="$(awk -F"\t" 'NR==3{print NF;exit}' ${INFILE})"

   if [ ${N_FIELDS} -gt 5 ];
   then
      awk '{print $1,$2,$3,$4,$7}' ${INFILE} > temp.txt
      cp -f temp.txt ${INFILE}
      rm temp.txt
   fi
    
if [ -f $OUTFILE ]; then
   echo "File $OUTFILE exists."
   echo "Appending to $OUTFILE."

   #Get number of tracks in files
   NEW_TRACKS="$(awk 'NR == 2' ${INFILE} | awk '{print $4}')"
   OLD_TRACKS="$(awk 'NR == 2' ${OUTFILE} | awk '{print $4}')"
   TOTAL_TRACKS=$((${NEW_TRACKS}+${OLD_TRACKS}))
   #Append new tracks to file
   awk 'NR > 3' ${INFILE} | awk -v N=${OLD_TRACKS} '{print ($1+N),$2,$3,$4,$5}' >> ${OUTFILE}
   #update track total
   NEW_HEADER="$(awk 'NR == 2' ${OUTFILE} | awk -v N=${TOTAL_TRACKS} '{print $1,$2,$3,N,$5}')"
   sed -i '.bak' '2s/.*/'"${NEW_HEADER}"'/' ${OUTFILE}
   rm ${OUTFILE}.bak
   
else
   echo "File $OUTFILE does not exist."
   echo "Creating $OUTFILE."
   cat ${INFILE} >> ${OUTFILE}
fi