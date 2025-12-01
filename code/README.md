# Code

The code in this folder should be the current working versions. 

When you submit a paper, commit the repo, and make a tag on that commit to say this was the version submitted to [journal] on [date]. 

We may add subdirectories for different steps of the analyses (but please do it while working, if not you lose relative paths). NB. don’t start script/folder names with numbers, use e.g. `s01_process_data.Rmd` instead.

R scripts and Rmarkdown will go here, and then they will write into the results folder. It is good for each script to write into a specific subdirectory, so it is easier to keep everything clean and understand which script created which output. 

Consider adding a command at the start of each script file to delete whatever is in the results folder for that script before re-making the results (can comment out if you don’t want to do this), e.g. r function `unlink()`.

Add a metafile that tells what every script does (e.g. you can add this info here in the readme).
