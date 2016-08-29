# start container and mount root
system("~/git/vol2bird/docker/start_container.sh /")
# pick a radarfile
radarfile="~/git/vol2bird/data/OPERA.h5"
# process it
data=system(paste("~/git/vol2bird/docker/run_vol2bird.sh",substring(normalizePath(radarfile),2)),intern=T)
# print the result
data
# stop the container
status=system("~/git/vol2bird/docker/stop_container.sh", ignore.stderr = T,intern=T)
# check that the container is stopped and removed
if(is.null(attributes(status)$status)) message("container stopped correctly...") else message("already no container running...")

