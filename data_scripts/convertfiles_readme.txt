#TR 20230531 Convert Files readme. This code was developed with the help of ChatGPT
#This script should be run from the same folder as process_date1.sh,etc. and gridfile.lcc
#Make sure you correctly point to the directories you want to access, too - the default is data/[directory]
#Then, run these commands if you haven't before in the session:
module load nco cdo
cd GEMMACH_data (or wherever you have all the folders and files)
#You only need to run these once per cluster, I believe (maybe rerun if you load a new version)
chmod +x process_date1.sh
chmod +x process_date2.sh
chmod +x process_date3.sh
#Sometimes getting out-of-memory faults on ncap2 - running with 10G/cpu works right now. 
#For setting your requirements, look at the node characteristics on the alliance website
#(e.g. https://docs.alliancecan.ca/wiki/Narval/en), and try to keep to one node. 
#Not sure how much this actually matters, but seems polite. For prioritization, try to use characteristics for more 
#common nodes (e.g. on Narval there are 1145 nodes with 250 GB of memory, but only 33 with 2009G)
#That being said, the files aren't requesting a full node so let her rip if you like