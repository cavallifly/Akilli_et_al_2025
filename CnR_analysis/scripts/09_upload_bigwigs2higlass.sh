echo "Consider all the analysed samples"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., cnr_H3K27me3_WD_NOnub_Rep1_L1_dm6_NA


if [[ $1 == "" ]];
then
    echo "Please provide a projectName to save your bigWigs in a higlass folder."
    exit
fi
projectName=$1

if [[ $2 != "" ]];
then
    replicates=$2
    echo "Provided sample to upload:"
    echo "${replicates}"
else
    replicates=$(ls -1 | grep cnr)
    echo ${replicates}
    echo "Please, provide the name of a file with the list of samples to upload."
    echo "Default: I am uploading all the .bigwigs in here ${PWD}"
fi

for replicate in ${replicates}
do
    echo ${replicate}
    
    assayName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $1}')
    targetName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $2}')
    cellName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $3}')
    condition=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $4}')
    replicateName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $5}')
    laneName=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $6}')
    assembly=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $7}')
    tag=$(echo ${replicate} | sed "s,_, ,g" | awk '{print $8}')
    
    cd ${replicate}

    for bwFile in $(ls -1 *RPKM*bw);
    do
	name=${bwFile%.bw}
	name=${name%.bigWig}
	name=${name%_normRPKM}	
	echo $bwFile $name

	docker exec higlass-container ls /tmp &> /dev/null

	check=$(docker exec higlass-container python higlass-server/manage.py list_tilesets | grep $name)
	if [[ $check != "" ]];
	then
	    echo "${bwFile} already available in the higlass"
	    echo "${check}"
	    echo "Remove the current entry using:"
	    echo "docker exec higlass-container python higlass-server/manage.py delete_tileset --uuid ${name}"
	    echo "If you wish to update it"
	    continue
	fi
	
	rsync -avz ${bwFile} /zdata/data/hg-tmp

	docker exec higlass-container python higlass-server/manage.py ingest_tileset --filename /tmp/${bwFile} --filetype bigwig --datatype vector --coordSystem ${assembly} --name "${name}" --project-name ${projectName} --uid "${name}"

	rm -v /zdata/data/hg-tmp/${bwFile}
    done

    cd ../ # Exit ${replicate}
    
    # Make for docker the tmp directory visible:
    #docker exec higlass-container ls /tmp
    # List tilesets:
    #docker exec higlass-container python higlass-server/manage.py list_tilesets
    # Delete tileset via uuid:
    #uuids=""
    #for uuid in ${uuids} ;
    #do
    #    docker exec higlass-container python higlass-server/manage.py delete_tileset --uuid ${uuid}
    #done
done
