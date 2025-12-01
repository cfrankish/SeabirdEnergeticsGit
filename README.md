# Template_project

This is the template project structure to help with project managment.  

This repository is set up as a `Template repository` which lets users generate new repositories with the same directory structure and files. 

Within the EvolEcolGroup Git, if you click `New` to create a new project repository you will be offered the chance to use a `Repository template; Start your repository with a template repository's contents`. 

Here there is a dropdown box where you can selcet this template `Template_project`. Download and personalise as needed.

A full tutorial is available at **update link**

### Back up strategy:

On this landing page detail your backup strategy for the project and provide links or 
comments on where these backups can be found.

## Basic rules to work on github

**Naming files and folders**
Never have spaces in your directory or file names (underscore can be used).
Please, only use lowercase letters when naming files and folders. 
Avoid naming your files starting with a number (e.g. instead of “01_process_data.R” use “s01_process_data.R”)
Do not use generic names such as “Figure_01”, “Figure_02”, etc. for your files and folders. Use more descriptive names as numbers change. It will make your life easier when you will need to work again on the project.
Naming your results folder as the script that created it (easy to find it after)  

**Paths**
You need to use absolute paths. 

**File size**
Small data (<1MB or so) can be kept on GitHub. For all the other data/files please use an alternate source (Google Drive, OneDrive, DropBox, etc.). We will work on scripts to retrieve data from different sources.  

**Integrating documents**
Documents are not properly read on GitHub. Create your document on google drive or a similar platform, and then link it within your README or any other markdown you are using.  

**Losing data**
GitHub does not delete anything; you can normally retrieve everything. 

**Commit**
Think carefully about your commit messages and branch names, as they will be very useful for returning to past changes for both you and others (for example when a project is made publicly available). Please make them informative.  

**Merge**
Once you merge a branch into master, kill the branch to avoid problems in the future.
Do not name a branch with the name of a branch that has been merged/deleted recently, as this may create problems. 

**Binary files**
Do not upload binary data (e.g.*.rds): changes in such files are too heavy to be handled by GitHub.  

## Github setup
If you already have a ssh key, follow these steps:
https://docs.github.com/en/authentication/connecting-to-github-with-ssh/checking-for-existing-ssh-keys 

If you need to generate a new key and add it to GitHub:
**1. Generate ssh key**
Source: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

**2. add key to account**
Source: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account  

**3. Test connection**
Source: https://docs.github.com/cn/authentication/connecting-to-github-with-ssh/testing-your-ssh-connection  

### Worked example (branch)  
**Summary**

