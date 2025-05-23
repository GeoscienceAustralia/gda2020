#!/bin/bash
#


#++++++++ get time for log file ++++++++
utc_yr4=`date -u +%Y`
#utc_yr4=2017
utc_month=`date -u +%m`
#utc_month=11
utc_mn_day=`date -u +%d`
#utc_mn_day=8
utc_hr_2=`date -u +%H`
#utc_hr_2=17

# minute
utc_min=`date -u +%M`

folder=$1
echo ${folder}

cd /data/craig
if [ ! -e "tmp" ]; then
   mkdir "tmp"
fi
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
metadata_file="RinexAntLs.txt"
usr_email="GDA2020@ga.gov.au"
#usr_email="alistair.deane@ga.gov.au"
#usr_email="jack.mccubbine@ga.gov.au"
#usr_email="Carl.Wang@ga.gov.au"

#++++++++ download antenne list from new AUSPOS web ++++++++
cd /data/craig

if [ ! -e "auspos_antenna" ]; then
   mkdir "auspos_antenna"
fi
cd /data/craig/auspos_antenna

#aws s3 cp s3://www.gaauspos.org/assets/data/antennaTypes.json ./
aws s3 cp s3://auspos-test-prod/assets/data/antennaTypes.json ./


#++++++++ go to folder in /data/craig/${folder}/rinexantls and get job list for the LOOP processing ++++++++
cd /data/craig/${folder}/rinexantls
rm /data/craig/${folder}/job_list.txt
ls -1 >> /data/craig/${folder}/job_list.txt 

input_list=/data/craig/${folder}/job_list.txt
while read line
do
    old_job=`echo -e "${line}" | awk '{print $1}'`
    echo "${old_job}"
    
    cd /data/craig/${folder}
    mkdir ${old_job}
    
    cp /data/craig/${folder}/rinexantls/${old_job} /data/craig/${folder}/${old_job}/${old_job}
    cp /data/craig/${folder}/rinexantls/${old_job} /data/craig/${folder}/${old_job}/${metadata_file}
    mv /data/craig/${folder}/rinexantls/${old_job} /data/craig/tmp/${metadata_file}
    
    no_rinex_file=`awk 'END {print NR}' /data/craig/tmp/${metadata_file}`
    for (( c=1; c<=${no_rinex_file}; c++ ))
    do
       usr_filename=`awk '{if(NR==ii) print $1}' ii="${c}" /data/craig/tmp/${metadata_file}`
       cp /data/craig/${folder}/${usr_filename} /data/craig/${folder}/${old_job}/${usr_filename}
       mv /data/craig/${folder}/${usr_filename} /data/craig/tmp/${usr_filename}
    done

    ##++++++++ download antenne list from new AUSPOS web ++++++++
    #cd /data/craig
    #
    #if [ ! -e "auspos_antenna" ]; then
    #   mkdir "auspos_antenna"
    #fi
    #cd /data/craig/auspos_antenna
    #
    #aws s3 cp s3://www.gaauspos.org/assets/data/antennaTypes.json ./
    #

    #++++++++ get last Craig's job number (string) from new AUSPOS web (job number start from 0000000000001) ++++++++
    last_job_str=`aws s3 ls s3://auspos-test-prod/uploads/  | sort | awk '{if(substr($2,14,23) == str_1) print substr($2,1,13)}' str_1="-00000000/" |tail -1`
    echo "last job string:"  ${last_job_str}
    
    #++++++++ crete new Craig's job number (number) from new AUSPOS web ++++++++ and transfer from number to string 0000000000002
    new_start_job_num=`expr ${last_job_str} + 1`
    new_start_job_str=`printf "%013i" "${new_start_job_num}"`
    echo "new start job string:"  ${new_start_job_str}
    short_new_start_job_str=`echo ${new_start_job_str} | cut -c10-13`
    #++++++++ start to create jason file +++++++++++++++
    cd /data/craig
    if [ ! -e "job" ]; then
       mkdir "job"
    fi
    
    cd /data/craig
    if [ ! -e "SNX" ]; then
       mkdir "SNX"
    fi

    cd /data/craig/job
    if [ ! -e "${new_start_job_str}" ]; then
       mkdir "${new_start_job_str}"
    else
       rm -r ${new_start_job_str}
       mkdir "${new_start_job_str}"
    fi
    
    path_base=/data/craig/job/${new_start_job_str}
    auspos_antenna_file="antennaTypes.json"
    rm /data/craig/auspos_antenna/tmp.txt
    
    num_auspos=`awk 'END {print NR}' /data/craig/auspos_antenna/${auspos_antenna_file}`
    num_auspos_sta=2
    num_auspos_end=`expr ${num_auspos} - 1`
    num_auspos_end_tmp=`expr ${num_auspos} - 2`
    
    for (( c=${num_auspos_sta}; c<=${num_auspos_end}; c++ ))
    do
       ante_meta=`awk '{if(NR==ii) print substr($0,4,20)}' ii="${c}" /data/craig/auspos_antenna/${auspos_antenna_file}`
       #ante_meta=`awk '{if(NR==ii) print substr($0,4,20) length(substr($0,4,23))}' ii="${c}" /data/craig/auspos_antenna/${auspos_antenna_file}`
       #length_ante_meta=`expr length "${ante_meta}"`
       #echo "ante_meta:"  ${length_ante_meta} >> /data/craig/auspos_antenna/1.txt
       echo ${ante_meta} >> /data/craig/auspos_antenna/tmp.txt
    done
    #
    #echo "num_auspos:"  ${num_auspos}
    #echo "num_auspos_sta:"  ${num_auspos_sta}
    #echo "num_auspos_end:"  ${num_auspos_end}
    #echo "num_auspos_end_tmp:"  ${num_auspos_end_tmp}
    
    
    num_rinex_file=`awk 'END {print NR}' /data/craig/tmp/${metadata_file}`
    echo "how many of RINEX files in metadata file:"  ${num_rinex_file}
    
    ouput_json_file=${path_base}/auspos_job.json
    echo "{" > $ouput_json_file
    echo "    \"files\": [" >> $ouput_json_file
    
    for (( c=1; c<=${num_rinex_file}; c++ ))
    do
       ante_meta_usr_1=`awk '{if(NR==ii) print $2}' ii="${c}" /data/craig/tmp/${metadata_file}`
       #ante_meta_usr_2=`awk '{if(NR==ii) print $3}' ii="${c}" ${path_base}/${metadata_file}`
       #if [ ${ante_meta_usr_2} != "NONE"]
       #echo "ante_meta_usr:"  ${ante_meta_usr_1}
       ante_meta_usr=${ante_meta_usr_1}
       ante_meta_usr_3=${ante_meta_usr_1}" NONE"
       echo "ante_meta_usr_3:"  ${ante_meta_usr_3}
       
       ante_meta_usr_file=`awk '{if(NR==ii) print $1}' ii="${c}" /data/craig/tmp/${metadata_file}`
       ante_meta_usr_height=`awk '{if(NR==ii) print $3}' ii="${c}" /data/craig/tmp/${metadata_file}`
    
       
       ante_meta_num_1=`awk '{if($0==amu3) print NR}' amu3="${ante_meta_usr_3}" /data/craig/auspos_antenna/tmp.txt`
       if  [[ -z $ante_meta_num_1 ]]; then
           echo "ANTENNA: " ${ante_meta_usr_3} " can not be found"
    	   echo "STOP!!"
    	   num_rinex_file="error"
       else
           echo "ante_meta:"  ${ante_meta_num_1} " " ${ante_meta_usr_3} >>${path_base}/tmp.txt
           ante_meta_num_2=`expr ${ante_meta_num_1} + 1`
    	   
    	   ante_meta_js=`awk '{if(NR==ii) print substr($0,4,23)}' ii="${ante_meta_num_2}" /data/craig/auspos_antenna/${auspos_antenna_file}`
    	   ante_meta_js=${ante_meta_usr_3}
    	   
           echo "        {" >> $ouput_json_file
           #echo "        \"antenna_type\" : \"${ante_meta_js}" >> $ouput_json_file
           echo "        \"antenna_type\" : \"${ante_meta_js}\"," >> $ouput_json_file
           echo "        \"antenna_height\" : \"${ante_meta_usr_height}\"," >> $ouput_json_file
           echo "        \"name\" : \"${ante_meta_usr_file}\"" >> $ouput_json_file
    	   
    	   if [ $c -eq ${num_rinex_file} ]; then
               echo "        }" >> $ouput_json_file
           else
               echo "        }," >> $ouput_json_file
           fi
    	   
    	   mv /data/craig/tmp/${ante_meta_usr_file} ${path_base}
       fi
    
       #for (( b=1; b<=${num_auspos_end_tmp}; b++ ))
       #do
       #   ante_meta=`awk '{if(NR==ii && $0==amu3) print $1 " " $2}' ii="${b}" amu3="${ante_meta_usr_3}" /data/craig/auspos_antenna/tmp.txt`
       #   echo "ante_meta:"  ${ante_meta} >>${path_base}/tmp.txt
       #   #if [ ${ante_meta_usr} == ${ante_meta} ]
       #   #ante_meta=`awk '{if(NR==ii) print substr($0,4,20) length(substr($0,4,23))}' ii="${c}" /data/craig/auspos_antenna/${auspos_antenna_file}`
       #   #length_ante_meta=`expr length "${ante_meta}"`
       #   #echo "ante_meta:"  ${length_ante_meta} >> /data/craig/auspos_antenna/1.txt
       #done
    done
    
    if  [ "${num_rinex_file}" != "error" ]; then
        cp /data/craig/tmp/${metadata_file} ${path_base}
    	
        echo "    ]," >> $ouput_json_file
        echo "    \"email\": \"${usr_email}\"," >> $ouput_json_file
        echo "    \"job_number\": \"${new_start_job_str}-00000000\"">> $ouput_json_file
        echo "}" >> $ouput_json_file
        
        #echo "${new_start_job_str}-00000000  ${num_rinex_file} -- ${folder} -- ${old_job} --- ${utc_yr4}-${utc_month}-${utc_mn_day}_${utc_hr_2}:${utc_min} -- https://d1hjwnsodwzynv.cloudfront.net/uploads/${new_start_job_str}-00000000/${short_new_start_job_str}_${new_start_job_str}-00000000.tar.gz -- https://d1hjwnsodwzynv.cloudfront.net/uploads/error_net/${short_new_start_job_str}_${new_start_job_str}-00000000.tar.gz" >> /data/craig/log_4_submit.txt
        echo "${new_start_job_str}-00000000  ${num_rinex_file} -- ${folder} -- ${old_job} --- ${utc_yr4}-${utc_month}-${utc_mn_day}_${utc_hr_2}:${utc_min} -- aws s3 cp s3://auspos-test-prod/uploads/${new_start_job_str}-00000000/${short_new_start_job_str}_${new_start_job_str}-00000000.tar.gz ./ -- aws s3 cp s3://auspos-test-prod/uploads/error_net/${short_new_start_job_str}_${new_start_job_str}-00000000.tar.gz ./ --" >> /data/craig/log_4_submit.txt

        echo "aws s3 cp s3://auspos-test-prod/uploads/${new_start_job_str}-00000000/ /data/craig/SNX --recursive  --exclude \"*\" --include \"*-*.SNX\" " >> /data/craig/down_4_SNX.sh
        #rm ${path_base}/tmp.txt
    
    	#for (( c=1; c<=${num_rinex_file}; c++ ))
        #do
    	#   file_name=`awk '{if(NR==ii) print $1}' ii="${c}" ${path_base}/${metadata_file}`
    	#done
    	mv /data/craig/tmp/${metadata_file} ${path_base}
        cp ${path_base}/auspos_job.json /data/craig/${folder}/${old_job}/
	
    	aws_s3_long_job=${new_start_job_str}-00000000
    	aws s3 cp ${path_base}/ s3://auspos-test-prod/uploads/${aws_s3_long_job}/  --recursive
    	
    	rm -f /data/craig/tmp/*

        sleep 2s

    fi

done < "${input_list}"
