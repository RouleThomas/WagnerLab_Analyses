# WagnerLab_Analyses

## Some usefull command ##

### Run command in background ###

```nohup [command] > log_files/log.txt 2>&1 &```
- ```&``` at the end of command will run it in the background
- log.files to create a log otherwise will appear in terminal once done
- ```nohup``` to maintain the command active even if you leave the server

```jobs -l``` display the status of all stopped and background jobs\
```fg``` brings a process to the foreground from the background\
```kill [process_id]``` kills the job in the background\



