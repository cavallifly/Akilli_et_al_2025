echo "Consider all the analysed replicates"
# Sample name assayName_targetName_cellName_condition_replicateName_laneName_assembly_tag
# e.g., rnaseq_RNA_WD_NOnub_Rep1_L1_dm6_NA
replicates=$(ls -1 | grep rnaseq | grep all)
echo "Checking ${replicates}"

assembly=$1
    
projectName=$2

#make for docker the tmp directory visible: 
docker exec higlass-container ls /tmp &> /dev/null

for replicate in ${replicates}
do
    echo ${replicate}

    cd $replicate
    files=$(ls -1 *bw *bigWig)
    #echo $files

    # List tilesets in higlass
    docker exec higlass-container python higlass-server/manage.py list_tilesets > _tmp_tilesets
    
    for file in $files ;
    do
	echo $file
	name=${file%.bw}
	name=${name%.bigWig}
	echo "Trying to load $file ($name) in higlass in the higlass project ${projectName}"

	checkTile=$(grep ${name} _tmp_tilesets | wc -l | awk '{print $1}')
	if [[ ${checkTile} -gt 0 ]];
	then
	    echo "The requested dataset ${name} is already present in the higlass visualizer"
	    grep ${name} _tmp_tilesets
	    continue
	fi

	rsync -avz ${file} /zdata/data/hg-tmp	
	docker exec higlass-container python higlass-server/manage.py ingest_tileset --filename /tmp/${file} --filetype bigwig --datatype vector --coordSystem ${assembly} --name "${name}" --project-name ${projectName} --uid "${name}"
	rm -v /zdata/data/hg-tmp/${file}
    done # Close cycle over $files
    cd ..
done # Close cycle over $replicates
exit
# Make for docker the tmp directory visible:
#docker exec higlass-container ls /tmp

# Delete tileset via uuid:
#uuids=""
#for uuid in ${uuids} ;
#do
#    docker exec higlass-container python higlass-server/manage.py delete_tileset --uuid ${uuid}
#done
