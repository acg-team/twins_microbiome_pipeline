# R wrapper for cut_adupters

# TODO : set path with command line argument

# run bash script
# https://stackoverflow.com/questions/11395217/run-a-bash-script-from-an-r-script
# form a bash string (make sure it has execution permissions)
x1 <- './0_cut_adapters.sh'   # script to run
arg <- ''  # add argument here (for the future)
x <- paste(x1, arg)
cat(x)  # sanity check

# run bash command
system(x)